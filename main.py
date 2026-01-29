import flet as ft
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from io import BytesIO
import base64
import threading
import os
import gdown

# ======================================================
# DATA MANAGER
# ======================================================
class DataManager:
    adata = None
    path = None
    pending_figs = []   # list of (fig, filename_suffix, fmt, dpi)

data = DataManager()

# ======================================================
# UTILITIES
# ======================================================
def fig_to_base64_str(fig, dpi=100):
    buf = BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=dpi)
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("utf-8")
    plt.close(fig)
    return b64

def save_fig(fig, path, fmt, dpi):
    try:
        fig.savefig(path, format=fmt, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
    except Exception as e:
        print(f"Error saving fig: {e}")

def glass(content, border_color=None, padding=20):
    return ft.Container(
        content=content,
        bgcolor=ft.Colors.WHITE,
        padding=padding,
        border_radius=18,
        border=ft.border.all(1, border_color) if border_color else None,
        shadow=ft.BoxShadow(
            blur_radius=15,
            spread_radius=1,
            color=ft.Colors.with_opacity(0.15, ft.Colors.BLACK),
        ),
    )

# ======================================================
# MAIN APP
# ======================================================
def main(page: ft.Page):

    PRIMARY = "#403EB2"
    SECONDARY = "#E91E63"
    BG = "#D4E7F6"
    SUCCESS_GREEN = "#4CAF50"
    ORANGE_LOADER = ft.Colors.ORANGE

    page.title = "scMalaria"
    page.window_width = 1300
    page.window_height = 850
    page.bgcolor = BG
    page.padding = 0

    # ======================================================
    # FILE PICKERS & SAVE LOGIC
    # ======================================================
    open_picker = ft.FilePicker()
    save_picker = ft.FilePicker()
    page.overlay.extend([open_picker, save_picker])

    upload_status_text = ft.Text("No file selected")
    project_info = ft.Text("-")

    # --- Local File Logic ---
    def on_file_picked(e):
        if e.files:
            data.path = e.files[0].path
            upload_status_text.value = f"Selected Local: {data.path}"
            load_local_btn.disabled = False
            page.update()

    open_picker.on_result = on_file_picked

    def on_save(e):
        if not e.path:
            return
        
        base_path = os.path.splitext(e.path)[0]
        count = 0
        for fig_data in data.pending_figs:
            fig = fig_data[0]
            name_suffix = fig_data[1]
            fmt = fig_data[2]
            dpi = fig_data[3]
            
            final_path = f"{base_path}_{name_suffix}.{fmt}"
            save_fig(fig, final_path, fmt, dpi)
            count += 1
        
        data.pending_figs.clear()
        page.snack_bar = ft.SnackBar(ft.Text(f"Saved {count} file(s) successfully!"))
        page.snack_bar.open = True
        page.update()

    save_picker.on_result = on_save

    # ======================================================
    # HELPER: DATATABLE GENERATOR
    # ======================================================
    def create_datatable(df: pd.DataFrame, title: str):
        # We slice head(50) for performance in rendering
        preview_df = df.head(50).reset_index()
        
        # Create Columns
        cols = []
        for col_name in preview_df.columns:
            cols.append(
                ft.DataColumn(
                    ft.Text(str(col_name), weight="bold", color="white"),
                )
            )

        # Create Rows
        rows = []
        for _, row in preview_df.iterrows():
            cells = []
            for idx, val in enumerate(row):
                # First column data styled (Index)
                if idx == 0:
                    cells.append(ft.DataCell(ft.Text(str(val), weight="bold", color=PRIMARY)))
                else:
                    cells.append(ft.DataCell(ft.Text(str(val), size=12, color="black")))
            rows.append(ft.DataRow(cells=cells))

        dt = ft.DataTable(
            columns=cols,
            rows=rows,
            heading_row_color=PRIMARY,
            border=ft.border.all(1, "grey"),
            vertical_lines=ft.border.BorderSide(1, "grey"),
            horizontal_lines=ft.border.BorderSide(1, "grey"),
        )

        return glass(
            ft.Column([
                ft.Text(f"{title} (Preview top 50 rows)", size=18, weight="bold", color=PRIMARY),
                ft.Container(
                    content=ft.Column([
                        ft.Row([dt], scroll=ft.ScrollMode.ALWAYS) # Horizontal Scroll
                    ], scroll=ft.ScrollMode.ALWAYS, height=300), # Vertical Scroll box
                    border=ft.border.all(1, "#eeeeee"),
                    border_radius=10,
                )
            ])
        )

    # Containers for dynamic upload content
    loading_container = ft.Column(
        visible=False,
        horizontal_alignment="center",
        controls=[
            ft.CupertinoActivityIndicator(radius=20, color=ORANGE_LOADER, animating=True),
            ft.Text("Downloading & Loading Data...", color=ORANGE_LOADER, weight="bold")
        ]
    )

    success_container = ft.Container(
        visible=False,
        bgcolor=SUCCESS_GREEN,
        padding=10,
        border_radius=8,
        content=ft.Row([
            ft.Icon(ft.Icons.CHECK_CIRCLE, color="white"),
            ft.Text("Uploaded Successful", color="white", weight="bold")
        ], alignment="center", spacing=10)
    )

    tables_container = ft.Column(visible=False, spacing=20)

    # ======================================================
    # HELPER: SAFE READ
    # ======================================================
    def safe_read_h5ad(filepath):
        """Reads h5ad and handles the 'gene_id' conflict error."""
        try:
            return sc.read_h5ad(filepath)
        except ValueError as ve:
            if "cannot insert gene_id, already exists" in str(ve):
                # Fallback: Read in 'r' mode then copy to memory to bypass index checks
                ad = sc.read_h5ad(filepath, backed='r')
                ad = ad.to_memory()
                return ad
            else:
                raise ve

    # ======================================================
    # LOAD DATA LOGIC
    # ======================================================
    def load_local_data():
        try:
            # UI Updates start
            loading_container.visible = True
            success_container.visible = False
            tables_container.visible = False
            upload_status_text.value = "Loading local file..."
            page.update()
            
            # Logic
            data.adata = safe_read_h5ad(data.path)
            
            # Sanity check for duplicate gene_id column
            if 'gene_id' in data.adata.var.columns and data.adata.var.index.name == 'gene_id':
                data.adata.var.rename(columns={'gene_id': 'gene_id_col'}, inplace=True)

            project_info.value = f"{data.adata.n_obs} cells × {data.adata.n_vars} genes"
            upload_status_text.value = "Success! Data loaded from Local Disk."
            
            # UI Updates End
            loading_container.visible = False
            success_container.visible = True
            
            # Build Tables
            tables_container.controls = [
                create_datatable(data.adata.obs, "Observation Metadata (adata.obs)"),
                create_datatable(data.adata.var, "Gene Metadata (adata.var)")
            ]
            tables_container.visible = True

        except Exception as ex:
            loading_container.visible = False
            upload_status_text.value = f"Error: {ex}"
        
        page.update()

    def load_drive_data(url):
        try:
            # UI Updates start
            loading_container.visible = True
            success_container.visible = False
            tables_container.visible = False
            upload_status_text.value = "Initiating Download..."
            page.update()
            
            output_file = "temp_drive_data.h5ad"
            
            # Download
            gdown.download(url, output_file, quiet=False, fuzzy=True)
            
            if not os.path.exists(output_file):
                loading_container.visible = False
                upload_status_text.value = "Error: Download failed or file not found."
                page.update()
                return

            upload_status_text.value = "Reading .h5ad file..."
            page.update()
            
            # Read
            data.adata = safe_read_h5ad(output_file)
            
            # Sanity check
            if 'gene_id' in data.adata.var.columns and data.adata.var.index.name == 'gene_id':
                data.adata.var.rename(columns={'gene_id': 'gene_id_col'}, inplace=True)

            project_info.value = f"{data.adata.n_obs} cells × {data.adata.n_vars} genes"
            upload_status_text.value = "Success! Data loaded from Cloud."
            
            # UI Updates End
            loading_container.visible = False
            success_container.visible = True
            
            # Build Tables
            tables_container.controls = [
                create_datatable(data.adata.obs, "Observation Metadata (adata.obs)"),
                create_datatable(data.adata.var, "Gene Metadata (adata.var)")
            ]
            tables_container.visible = True
            
        except Exception as ex:
            loading_container.visible = False
            upload_status_text.value = f"Drive Error: {ex}"
            print(f"Drive Error Details: {ex}")
        
        page.update()

    # --- Upload View Controls ---
    
    # 1. Local
    browse_btn = ft.ElevatedButton(
        "Browse Local .h5ad",
        icon=ft.Icons.FOLDER_OPEN,
        on_click=lambda e: open_picker.pick_files(allowed_extensions=["h5ad"]),
    )
    load_local_btn = ft.ElevatedButton(
        "Load Local",
        icon=ft.Icons.UPLOAD_FILE,
        disabled=True,
        on_click=lambda e: threading.Thread(target=load_local_data, daemon=True).start(),
    )

    # 2. Cloud (Drive)
    default_url = "https://drive.google.com/file/d/1EyMVVP5FOUa7jzLnMDcIerZEEIprHu3m/view?usp=sharing"
    drive_url_field = ft.TextField(
        label="Google Drive URL", 
        value=default_url, 
        text_size=12,
        expand=True
    )
    
    load_drive_btn = ft.ElevatedButton(
        "Load",
        icon=ft.Icons.CLOUD_DOWNLOAD,
        bgcolor=SECONDARY,
        color="white",
        on_click=lambda e: threading.Thread(target=load_drive_data, args=(drive_url_field.value,), daemon=True).start()
    )

    upload_view = ft.Column(
        [
            glass(
                ft.Column(
                    [
                        ft.Row([ft.Icon(ft.Icons.CLOUD_UPLOAD, color=PRIMARY, size=30), ft.Text("Data Source", size=24, weight="bold", color=PRIMARY)]),
                        ft.Divider(),
                        
                        # Cloud Section
                        ft.Text("Option 1: Cloud (Default)", weight="bold", color=SECONDARY),
                        ft.Row([drive_url_field, load_drive_btn]),
                        
                        # Loading & Success Indicators
                        loading_container,
                        success_container,

                        ft.Divider(),
                        
                        # Local Section
                        ft.Text("Option 2: Local Disk", weight="bold"),
                        ft.Row([browse_btn, load_local_btn]),
                        
                        ft.Divider(),
                        ft.Text("Log:", weight="bold", size=12),
                        upload_status_text,
                    ],
                    spacing=15,
                )
            ),
            # Table Section
            tables_container
        ],
        scroll=ft.ScrollMode.AUTO,
        expand=True,
        spacing=20
    )

    # ======================================================
    # EXPLORE VIEW
    # ======================================================
    def explore_view():
        if data.adata is None:
            return glass(ft.Text("Please upload data first"))

        adata = data.adata

        # --- Metadata Logic ---
        meta_plot = ft.Image(visible=False, width=600, height=450, fit=ft.ImageFit.CONTAIN)
        meta_dd = ft.Dropdown(
            label="Metadata",
            options=[ft.dropdown.Option(c) for c in adata.obs.columns if adata.obs[c].nunique() <= 50],
        )
        meta_fmt = ft.Dropdown(label="Format", value="png", options=[ft.dropdown.Option(f) for f in ["png", "pdf", "svg", "jpeg"]], width=80)
        meta_w = ft.TextField(label="W", value="6", width=60)
        meta_h = ft.TextField(label="H", value="6", width=60)
        meta_dpi = ft.TextField(label="DPI", value="300", width=70)

        def plot_meta(e):
            def task():
                try:
                    fig, ax = plt.subplots(figsize=(6, 5), dpi=100)
                    sc.pl.umap(adata, color=meta_dd.value, show=False, ax=ax, frameon=False, title=meta_dd.value)
                    meta_plot.src_base64 = fig_to_base64_str(fig)
                    meta_plot.visible = True
                    page.update()
                except Exception as ex:
                    print(ex)
            threading.Thread(target=task, daemon=True).start()

        def download_meta(e):
            try:
                w = float(meta_w.value)
                h = float(meta_h.value)
                d = int(meta_dpi.value)
                fmt = meta_fmt.value
                
                fig, ax = plt.subplots(figsize=(w, h), dpi=d)
                sc.pl.umap(adata, color=meta_dd.value, show=False, ax=ax, frameon=False, title=meta_dd.value)
                
                data.pending_figs = [(fig, f"meta_{meta_dd.value}", fmt, d)]
                save_picker.save_file(file_name=f"umap_{meta_dd.value}.{fmt}")
            except Exception as ex:
                page.snack_bar = ft.SnackBar(ft.Text(f"Error preparing download: {ex}"))
                page.snack_bar.open = True
                page.update()

        meta_section = glass(
            ft.Column([
                ft.Text("Metadata UMAP", size=20, weight="bold", color=PRIMARY),
                ft.Row([
                    ft.Column([meta_dd, ft.ElevatedButton("Update Plot", bgcolor=PRIMARY, color="white", on_click=plot_meta)], width=250),
                    meta_plot
                ], alignment=ft.MainAxisAlignment.START, vertical_alignment=ft.CrossAxisAlignment.START),
                ft.Divider(),
                ft.Row([
                    ft.Text("Export:", weight="bold"),
                    meta_w, meta_h, meta_dpi, meta_fmt, 
                    ft.IconButton(ft.Icons.DOWNLOAD, on_click=download_meta, tooltip="Download with these settings")
                ], vertical_alignment=ft.CrossAxisAlignment.CENTER)
            ])
        )

        # --- Gene Logic ---
        genes = sorted(adata.var_names.tolist())
        if "symbol_merged" in adata.var.columns:
             genes = sorted(adata.var["symbol_merged"].dropna().unique().tolist())
        
        selected_genes = set()
        plot_grid = ft.Column(spacing=20)
        gene_search = ft.TextField(label="Search gene", prefix_icon=ft.Icons.SEARCH)
        gene_checklist = ft.Column(height=200, scroll=ft.ScrollMode.AUTO)
        
        gene_fmt = ft.Dropdown(label="Format", value="png", options=[ft.dropdown.Option(f) for f in ["png", "pdf", "svg", "jpeg"]], width=80)
        gene_w = ft.TextField(label="W", value="5", width=60)
        gene_h = ft.TextField(label="H", value="5", width=60)
        gene_dpi = ft.TextField(label="DPI", value="300", width=70)

        def refresh_gene_list():
            gene_checklist.controls.clear()
            q = gene_search.value.lower()
            count = 0
            for g in genes:
                if q in g.lower():
                    gene_checklist.controls.append(
                        ft.Checkbox(label=g, value=g in selected_genes, on_change=lambda e, gene=g: toggle_gene(gene, e.control.value))
                    )
                    count += 1
                    if count > 100: break 
            page.update()

        def toggle_gene(gene, checked):
            if checked:
                if len(selected_genes) >= 4:
                    page.snack_bar = ft.SnackBar(ft.Text("Max 4 genes allowed"))
                    page.snack_bar.open = True
                    page.update()
                    refresh_gene_list()
                    return
                selected_genes.add(gene)
            else:
                selected_genes.discard(gene)

        gene_search.on_change = lambda e: refresh_gene_list()
        layer_dd = ft.Dropdown(label="Layer", value="log1p" if "log1p" in adata.layers else None, options=[ft.dropdown.Option(l) for l in adata.layers.keys()])
        cmap_dd = ft.Dropdown(label="Cmap", value="magma", options=[ft.dropdown.Option(c) for c in ["viridis", "plasma", "magma", "RdBu_r"]])

        def plot_genes(e):
            def task():
                plot_grid.controls.clear()
                blocks = []
                for gene in selected_genes:
                    try:
                        fig, ax = plt.subplots(figsize=(5, 5), dpi=100)
                        if "symbol_merged" in adata.var.columns:
                            sc.pl.umap(adata, color=gene, gene_symbols="symbol_merged", layer=layer_dd.value, cmap=cmap_dd.value, show=False, ax=ax, frameon=False)
                        else:
                            sc.pl.umap(adata, color=gene, layer=layer_dd.value, cmap=cmap_dd.value, show=False, ax=ax, frameon=False)
                        
                        blocks.append(ft.Column([ft.Text(gene, weight="bold", color=PRIMARY), ft.Image(src_base64=fig_to_base64_str(fig), width=300, height=300)], horizontal_alignment="center"))
                    except Exception as err:
                        print(f"Error plotting {gene}: {err}")

                for i in range(0, len(blocks), 2):
                    plot_grid.controls.append(ft.Row(blocks[i:i + 2], spacing=30))
                page.update()

            threading.Thread(target=task, daemon=True).start()

        def download_genes(e):
            if not selected_genes: return
            w = float(gene_w.value)
            h = float(gene_h.value)
            d = int(gene_dpi.value)
            fmt = gene_fmt.value
            data.pending_figs = []
            
            for gene in selected_genes:
                fig, ax = plt.subplots(figsize=(w, h), dpi=d)
                if "symbol_merged" in adata.var.columns:
                    sc.pl.umap(adata, color=gene, gene_symbols="symbol_merged", layer=layer_dd.value, cmap=cmap_dd.value, show=False, ax=ax, frameon=False)
                else:
                    sc.pl.umap(adata, color=gene, layer=layer_dd.value, cmap=cmap_dd.value, show=False, ax=ax, frameon=False)
                
                data.pending_figs.append((fig, f"gene_{gene}", fmt, d))
            
            save_picker.save_file(file_name=f"genes_multi.{fmt}")

        gene_section = glass(
            ft.Column([
                ft.Text("Multi-Gene Expression", size=20, weight="bold", color=PRIMARY),
                ft.Row([
                    ft.Column([gene_search, gene_checklist, layer_dd, cmap_dd, ft.ElevatedButton("Plot Genes", bgcolor=PRIMARY, color="white", on_click=plot_genes)], width=300),
                    plot_grid
                ], vertical_alignment=ft.CrossAxisAlignment.START),
                ft.Divider(),
                ft.Row([
                    ft.Text("Export:", weight="bold"),
                    gene_w, gene_h, gene_dpi, gene_fmt,
                    ft.ElevatedButton("Download All", icon=ft.Icons.DOWNLOAD, on_click=download_genes)
                ])
            ])
        )
        
        refresh_gene_list()
        return ft.Column([meta_section, gene_section], spacing=20, scroll=ft.ScrollMode.AUTO)

    # ======================================================
    # FUNCTIONAL VIEW
    # ======================================================
    def functional_view():
        if data.adata is None:
            return glass(ft.Text("Please upload data first"))
        
        adata = data.adata
        fd = adata.uns.get('plasmoDB_functional_analysis_pb', {})
        
        if not fd:
            return glass(ft.Column([
                ft.Icon(ft.Icons.WARNING, color="orange", size=40),
                ft.Text("No functional analysis data found.", size=16)
            ]))

        current_df = [None]
        current_viz = [None]
        
        # --- Controls ---
        embs = [k for k in adata.obsm.keys() if "X_" in k]
        def_emb = "X_umap" if "X_umap" in embs else embs[0]
        dd_basis = ft.Dropdown(label="Basis", options=[ft.dropdown.Option(k) for k in embs], value=def_emb)
        
        times = sorted(list(fd.keys()))
        dd_time = ft.Dropdown(label="Timepoint", options=[ft.dropdown.Option(t) for t in times], value=times[0] if times else None)
        dd_cat = ft.Dropdown(label="Category")
        slider_top = ft.Slider(min=5, max=50, divisions=45, value=15, label="Top: {value}")
        dd_pathway = ft.Dropdown(label="Select Pathway", disabled=True)
        
        # --- Export Settings ---
        func_w = ft.TextField(label="W", value="6", width=60)
        func_h = ft.TextField(label="H", value="6", width=60)
        func_dpi = ft.TextField(label="DPI", value="300", width=70)
        func_fmt = ft.Dropdown(label="Fmt", value="png", options=[ft.dropdown.Option(f) for f in ["png", "pdf", "svg"]], width=80)
        
        export_expander = ft.ExpansionTile(
            title=ft.Text("Export Settings", size=14, color=SECONDARY),
            controls=[
                ft.Row([func_w, func_h, func_dpi, func_fmt], alignment=ft.MainAxisAlignment.CENTER)
            ]
        )

        summary_container = ft.Column()
        score_container = ft.Column()
        
        def update_cats(e=None):
            t = dd_time.value
            if t and t in fd:
                cats = sorted(list(fd[t].keys()))
                dd_cat.options = [ft.dropdown.Option(c) for c in cats]
                dd_cat.value = cats[0] if cats else None
            page.update()

        update_cats()
        dd_time.on_change = update_cats

        def generate_summary(e):
            t = dd_time.value
            c = dd_cat.value
            n_p = int(slider_top.value)
            
            if t and c:
                df = fd[t][c].copy()
                pval_col = next((col for col in df.columns if 'P-value' in col or 'Benjamini' in col), None)
                term_col = next((col for col in df.columns if 'Name' in col or 'Term' in col), None)
                
                if pval_col and term_col:
                    df[pval_col] = pd.to_numeric(df[pval_col], errors='coerce').fillna(1.0)
                    df['neg_log_p'] = -np.log10(df[pval_col].replace(0, 1e-300))
                    
                    viz = df.sort_values('neg_log_p', ascending=False).head(n_p).reset_index(drop=True)
                    current_df[0] = df
                    current_viz[0] = viz
                    
                    dd_pathway.options = [ft.dropdown.Option(row[term_col]) for _, row in viz.iterrows()]
                    dd_pathway.value = viz.iloc[0][term_col] if not viz.empty else None
                    dd_pathway.disabled = False
                    
                    fig, ax = plt.subplots(figsize=(8, 6))
                    y_pos = np.arange(len(viz))
                    ax.barh(y_pos, viz['neg_log_p'], color=plt.cm.viridis(viz['neg_log_p'] / viz['neg_log_p'].max()))
                    ax.set_yticks(y_pos)
                    ax.set_yticklabels(viz[term_col])
                    ax.invert_yaxis()
                    ax.set_xlabel('-log10(P-value)')
                    ax.set_title(f"Top {n_p} Pathways ({c} - {t})")
                    
                    summary_container.controls = [
                        glass(
                            ft.Column([
                                ft.Text("Pathway Enrichment", weight="bold", size=18, color=PRIMARY),
                                ft.Image(src_base64=fig_to_base64_str(fig), width=650, height=500, fit=ft.ImageFit.CONTAIN)
                            ], horizontal_alignment="center")
                        )
                    ]
                    page.update()

        def save_functional_umap(e):
            p_name = dd_pathway.value
            try:
                w = float(func_w.value)
                h = float(func_h.value)
                d = int(func_dpi.value)
                fmt = func_fmt.value
                
                fig_u, ax_u = plt.subplots(figsize=(w, h), dpi=d)
                sc.pl.embedding(adata, basis=dd_basis.value, color='pathway_score', cmap='RdBu_r', 
                                title=f"{p_name}", frameon=False, show=False, ax=ax_u)
                
                data.pending_figs = [(fig_u, f"score_{p_name[:15]}", fmt, d)]
                save_picker.save_file(file_name=f"pathway_score.{fmt}")
            except Exception as ex:
                print(ex)

        def compute_score(e):
            p_name = dd_pathway.value
            viz = current_viz[0]
            
            if p_name and viz is not None:
                term_col = next((col for col in viz.columns if 'Name' in col or 'Term' in col), None)
                genes_col = next((col for col in viz.columns if 'gene list' in col or 'Genes' in col), None)
                
                row = viz[viz[term_col] == p_name].iloc[0]
                raw_g = str(row[genes_col]).replace('"','').replace("'","").replace("[","").replace("]","")
                gl = [g.strip() for g in (raw_g.split(",") if "," in raw_g else raw_g.split())]
                v_gl = [g for g in gl if g in adata.var_names]
                
                if not v_gl:
                    page.snack_bar = ft.SnackBar(ft.Text("No valid genes found."))
                    page.snack_bar.open = True
                    page.update()
                    return

                score_container.controls = [ft.ProgressBar(width=200, color=SECONDARY)]
                page.update()
                
                def task():
                    try:
                        sc.tl.score_genes(adata, v_gl, score_name='pathway_score')
                        
                        fig_u, ax_u = plt.subplots(figsize=(5, 5))
                        sc.pl.embedding(adata, basis=dd_basis.value, color='pathway_score', cmap='RdBu_r', 
                                        title=f"{p_name}", frameon=False, show=False, ax=ax_u)
                        
                        labels = []
                        for g in v_gl:
                            if "symbol" in adata.var.columns: 
                                labels.append(adata.var.loc[g, "symbol"] if isinstance(adata.var.loc[g, "symbol"], str) else g)
                            else: 
                                labels.append(g)
                        
                        # Flatten logic
                        raw_mean = adata[:, v_gl].X.mean(axis=0)
                        mean_expr = np.array(raw_mean).flatten() 

                        dfe = pd.DataFrame({'Gene': labels, 'Expr': mean_expr}).sort_values('Expr', ascending=True).tail(15)
                        
                        fig_b, ax_b = plt.subplots(figsize=(5, 5))
                        ax_b.barh(dfe['Gene'], dfe['Expr'], color='#D8BFD8')
                        ax_b.set_xlabel("Mean Expression")
                        ax_b.set_title("Top Driver Genes")

                        score_container.controls = [
                            glass(
                                ft.Column([
                                    ft.Text("Pathway Activity & Drivers", weight="bold", size=18, color=PRIMARY),
                                    ft.Row([
                                        ft.Column([
                                            ft.Text("Activity UMAP", weight="bold"),
                                            ft.Image(src_base64=fig_to_base64_str(fig_u), width=400, height=400),
                                            ft.ElevatedButton("Save UMAP", icon=ft.Icons.SAVE, on_click=save_functional_umap)
                                        ], horizontal_alignment="center"),
                                        ft.Column([
                                            ft.Text("Top Genes", weight="bold"),
                                            ft.Image(src_base64=fig_to_base64_str(fig_b), width=400, height=400)
                                        ], horizontal_alignment="center")
                                    ], alignment=ft.MainAxisAlignment.CENTER, spacing=30)
                                ])
                            )
                        ]
                        
                    except Exception as ex:
                        score_container.controls = [ft.Text(f"Error: {ex}", color="red")]
                        print(ex) 
                    
                    page.update()
                
                threading.Thread(target=task, daemon=True).start()

        controls_panel = glass(
            ft.Column([
                ft.Text("Settings", size=18, weight="bold", color=SECONDARY),
                dd_basis,
                dd_time,
                dd_cat,
                slider_top,
                ft.ElevatedButton("Show Pathways", bgcolor=SECONDARY, color="white", on_click=generate_summary),
                ft.Divider(),
                dd_pathway,
                ft.ElevatedButton("Compute Score", bgcolor=ft.Colors.BLACK87, color="white", icon=ft.Icons.CALCULATE, on_click=compute_score),
                ft.Divider(),
                export_expander
            ], spacing=10),
            border_color=SECONDARY
        )

        return ft.Row(
            [
                ft.Column([controls_panel], width=300),
                ft.Column([summary_container, score_container], expand=True, scroll=ft.ScrollMode.AUTO)
            ],
            vertical_alignment=ft.CrossAxisAlignment.START,
            spacing=20
        )

    # ======================================================
    # TABS & LAYOUT
    # ======================================================
    sidebar = ft.Container(
        width=250, bgcolor=ft.Colors.WHITE, padding=20,
        border=ft.border.only(right=ft.border.BorderSide(1, "#E0E0E0")),
        content=ft.Column([
            ft.Text("Project Info", weight="bold", size=16, color="grey"),
            ft.Divider(),
            ft.Row([ft.Icon(ft.Icons.DATASET, size=16), project_info]),
            ft.Divider(),
            ft.Text("Status", weight="bold", size=12),
            upload_status_text
        ])
    )

    body = ft.Container(expand=True, padding=20)

    def tab_change(e):
        idx = e.control.selected_index
        if idx == 0: body.content = upload_view
        elif idx == 1: body.content = explore_view()
        elif idx == 2: body.content = functional_view()
        page.update()

    tabs = ft.Tabs(
        on_change=tab_change, label_color="white", unselected_label_color="white60",
        indicator_color=SECONDARY, divider_color="transparent",
        tabs=[
            ft.Tab(text="Upload", icon=ft.Icons.UPLOAD_FILE),
            ft.Tab(text="Explore", icon=ft.Icons.EXPLORE),
            ft.Tab(tab_content=ft.Row([
                ft.Icon(ft.Icons.SCIENCE), ft.Text("Functional (Liver-stage)"),
                ft.Container(content=ft.Text("NEW", size=9, color="white", weight="bold"), bgcolor=SECONDARY, padding=ft.padding.symmetric(horizontal=4), border_radius=4)
            ])),
        ],
    )

    topbar = ft.Container(
        gradient=ft.LinearGradient(colors=[PRIMARY, "#5E5CCE"]),
        padding=ft.padding.symmetric(horizontal=20, vertical=10),
        shadow=ft.BoxShadow(blur_radius=5, color=ft.Colors.BLACK26),
        content=ft.Column([
            ft.Row([ft.Icon(ft.Icons.BIOTECH_ROUNDED, color="white", size=32), ft.Text("scMalaria", size=24, weight="bold", color="white")]),
            tabs,
        ]),
    )

    page.add(ft.Column([topbar, ft.Row([sidebar, body], expand=True)], expand=True))
    body.content = upload_view
    page.update()

if __name__ == "__main__":
    ft.app(target=main)

    