import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from .id_normalizer import normalize_genome_id, format_genome_display
from .drep_cluster_func import prepare_cluster_data, simplify_pathway_class


def select_top_modules(df: pd.DataFrame, top_modules: int | None = None, mode: str = "mean") -> pd.DataFrame:
    """
    Select the most relevant pathway modules for plotting.
    This function ranks modules across all genomes using:
    - mean: modules with the highest average completeness (core functions)
    - variance: module with the largest variation in completeness (differential functions)
    
    Returns filtered dataframe containing only selected modules.
    
    """
    if top_modules is None:
        return df.copy()
    
    # Build matrix
    # rows = genomes, columns = modules, values = completeness
    pivot = df.pivot_table(
        index="contig",
        columns="module_accession",
        values="completeness",
        aggfunc="max",
    )

    pivot = pivot.fillna(0)

    if mode == "mean":
        scores = pivot.mean(axis=0)
    elif mode == "variance":
        scores = pivot.var(axis=0)
    else:
        raise ValueError("mode must be 'mean' or 'variance'")

    keep_modules = scores.sort_values(ascending=False).head(top_modules).index
    return df[df["module_accession"].isin(keep_modules)].copy()


def simplify_module_name(name: str) -> str:
    if pd.isna(name):
        return ""

    name = str(name)

    # cut at "("
    if "(" in name:
        name = name.split("(")[0]

    # cut at ","
    if "," in name:
        name = name.split(",")[0]
    return name.strip()


def build_module_labels(module_meta: pd.DataFrame, module_label: str = "id") -> list[str]:
    """
    x-axis labels for KEGG modules
    module_label:
    - id: module accession
    - name: module name
    - both: accession and name
    """
    if module_label == "id":
        return module_meta["module_accession"].astype(str).tolist()

    if module_label == "name":
        return module_meta["pathway_name"].apply(simplify_module_name).tolist()

    if module_label == "both":
        return (
            module_meta["module_accession"].astype(str)
            + " | "
            + module_meta["pathway_name"].apply(simplify_module_name)
        ).tolist()
    raise ValueError("module_label must be one of: id, name, both")


def prepare_heatmap_data(pathway_df: pd.DataFrame, top_modules: int | None = None,
    mode: str = "mean", representatives_df: pd.DataFrame = None,
    gtdb_df: pd.DataFrame = None, top_representatives: int | None = None):
    """
    Prepare all data needed to draw the pathway module heatmap.
    Build completeness matrix and pathway class block for heatmap.
    """
    df = pathway_df.copy()

    required_cols = {"contig", "module_accession", "completeness", "pathway_name", "pathway_class",}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"pathway file missing required columns: {sorted(missing)}")

    df["contig_normalized"] = df["contig"].apply(normalize_genome_id)
    df["completeness"] = pd.to_numeric(df["completeness"], errors="coerce")
    df["pathway_class_simple"] = df["pathway_class"].apply(simplify_pathway_class)

    representatives_only = top_representatives is not None

    if representatives_only:
        if representatives_df is None or gtdb_df is None:
            raise ValueError(
                "representatives_only=True requires both representatives_df (dRep) and gtdb_df"
            )
        
        pw_genomes = set(df["contig_normalized"].unique())

        # Build cluster data and keep larger set first
        cluster_data, _, _ = prepare_cluster_data(
            drep_df=representatives_df,
            gtdb_df=gtdb_df,
            tax_levels=("phylum",),
            top_n=top_representatives * 2,
            pathway_df=pathway_df,
        )

        cluster_data = cluster_data[
            cluster_data["representative"].isin(pw_genomes)
            ].reset_index(drop=True)

        # prepare_cluster_data sorts in ascending order, keep largest clusters with tail()
        cluster_data = cluster_data.tail(top_representatives).reset_index(drop=True)

        representatives = cluster_data["representative"].tolist()
        row_order = [r for r in representatives if r in set(df["contig_normalized"])]
        rep_to_size = dict(zip(cluster_data["representative"], cluster_data["n_members"]))
    else:
        # If all MAGs, sort them alphabetically
        row_order = sorted(df["contig_normalized"].dropna().unique().tolist())
        rep_to_size = {}

    df = select_top_modules(df, top_modules=top_modules, mode=mode)

    # Sort by pathway class
    module_meta = (
        df[["module_accession", "pathway_name", "pathway_class_simple"]]
        .drop_duplicates()
        .sort_values(["pathway_class_simple", "module_accession"])
        .reset_index(drop=True)
    )

    module_order = module_meta["module_accession"].tolist()

    # Build heatmap matrix: rows = genomes, columns = modules, values = completeness
    heatmap_df = (
        df.pivot_table(
            index="contig_normalized",
            columns="module_accession",
            values="completeness",
            aggfunc="max",
        )
        .reindex(index=row_order, columns=module_order)
    )

    # identify pathway-class blocks for grouping
    class_blocks = []
    if not module_meta.empty:
        start = 0
        current_class = module_meta.loc[0, "pathway_class_simple"]

        for i in range(1, len(module_meta)):
            cls = module_meta.loc[i, "pathway_class_simple"]
            if cls != current_class:
                class_blocks.append((current_class, start, i - 1))
                start = i
                current_class = cls
        class_blocks.append((current_class, start, len(module_meta) - 1))

    cluster_sizes = [rep_to_size.get(rep, np.nan) for rep in heatmap_df.index]

    return heatmap_df, module_meta, class_blocks, cluster_sizes


def pathway_module_heatmap(pathway_df: pd.DataFrame, output_path: str, top_modules: int | None = None,
    mode: str = "mean", fmt: str = "png", fig_size=None, row_fontsize: int = 8, 
    module_label: str = "id", representatives_df: pd.DataFrame = None, 
    gtdb_df: pd.DataFrame = None, top_representatives: int | None = None):
    """
    Create a heatmap shwoing pathway/module completeness across MAGs.

    - rows: MAGs / contigs
    - columns: module_accession
    - values: pathway completeness (0% - 100%)
    """
    os.makedirs(output_path, exist_ok=True)

    representatives_only = top_representatives is not None

    heatmap_df, module_meta, class_blocks, cluster_sizes = prepare_heatmap_data(
        pathway_df=pathway_df,
        top_modules=top_modules,
        mode=mode,
        representatives_df=representatives_df,
        gtdb_df=gtdb_df,
        top_representatives=top_representatives,
    )

    if heatmap_df.empty:
        raise ValueError("Heatmap matrix is empty")

    n_rows, n_cols = heatmap_df.shape

    if fig_size is None:
        width = max(12, 5 + n_cols * 0.28)
        height = max(8, 3 + n_rows * 0.22)
        fig_size = (width, height)

    fig, ax = plt.subplots(figsize=fig_size)

    # Flip row order so the largest cluster appears at the top
    heatmap_df = heatmap_df.iloc[::-1]
    cluster_sizes = cluster_sizes[::-1]

    values = heatmap_df.to_numpy(dtype=float)
    masked = np.ma.masked_invalid(values)

    cmap = plt.get_cmap("Reds").copy()
    cmap.set_bad(color="white")
    norm = Normalize(vmin=0, vmax=100)

    im = ax.imshow(masked, aspect="auto", interpolation="none", cmap=cmap, norm=norm)

    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)

    # MAG labels on the y-axis
    ax.set_yticks(np.arange(n_rows))
    display_names = [format_genome_display(x) for x in heatmap_df.index]
    ax.set_yticklabels(display_names, fontsize=row_fontsize)

    # Module labels on x-axis if wanted
    ax.set_xticks(np.arange(n_cols))

    module_labels = build_module_labels(module_meta, module_label=module_label)

    ax.set_xticklabels(
        module_labels,
        rotation=90,
        fontsize=8 if module_label in {"name", "both"} or n_cols > 30 else 8,
    )

    ax.tick_params(axis="x", bottom=True, labelbottom=True, top=False, length=2)

    # Grid lines between heatmap cells
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="#bdbdbd", linewidth=0.35)
    ax.tick_params(which="minor", bottom=False, left=False)

    for spine in ax.spines.values():
        spine.set_visible(False)

    # separator lines between pathway classes
    for _, _, end in class_blocks[:-1]:
        ax.axvline(end + 0.5, color="white", linewidth=2.2)

    # colored bars & staggered labels with connector lines
    class_palette = [
        "#4C78A8", "#F58518", "#54A24B", "#E45756", "#72B7B2",
        "#B279A2", "#FF9DA6", "#9D755D", "#BAB0AC", "#2E91E5"
    ]

    y_bar = -0.95
    bar_height = 0.22
    label_levels = [-1.55, -2.25, -2.95, -3.65]
    line_bottom_pad = 0.01

    # Colored bars and class labels above heatmap
    for i, (cls, start, end) in enumerate(class_blocks):
        center = (start + end) / 2.0
        width = (end - start) + 1
        color = class_palette[i % len(class_palette)]

        rect = plt.Rectangle(
            (start - 0.5, y_bar),
            width,
            bar_height,
            facecolor=color,
            edgecolor="white",
            linewidth=0.8,
            clip_on=False,
            zorder=4,
        )
        ax.add_patch(rect)

        y_text = label_levels[i % len(label_levels)]

        ax.plot(
            [center, center],
            [y_bar, y_text - line_bottom_pad],
            color=color,
            linewidth=1.0,
            clip_on=False,
            zorder=5,
        )

        ax.text(
            center,
            y_text,
            cls,
            ha="center",
            va="bottom",
            fontsize=9,
            color=color,
            clip_on=False,
            zorder=6,
        )
    
    # ---- vertical side heatmap: cluster size / MAG number ----
    size_im = None
    if representatives_only and len(cluster_sizes) == n_rows:
        pos = ax.get_position()

        size_vals = np.array(cluster_sizes, dtype=float).reshape(-1, 1)

        size_ax = fig.add_axes([pos.x1 + 0.01, pos.y0, 0.02, pos.height])

        size_cmap = plt.get_cmap("Oranges").copy()
        size_cmap.set_bad(color="white")

        size_im = size_ax.imshow(
            size_vals,
            aspect="auto",
            interpolation="none",
            cmap=size_cmap,
        )

        size_ax.set_ylim(n_rows - 0.5, -0.5)
        size_ax.set_xticks([])
        size_ax.set_yticks([])

        size_ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
        size_ax.grid(which="minor", color="#bdbdbd", linewidth=0.35)
        size_ax.tick_params(which="minor", left=False, bottom=False)

        for spine in size_ax.spines.values():
            spine.set_edgecolor("#8c8c8c")
            spine.set_linewidth(0.8)

        size_ax.text(
            0.5,
            1.02,
            "MAG\ncount",
            ha="center",
            va="bottom",
            fontsize=8,
            transform=size_ax.transAxes,
        )

    title_text = {
        "mean": "MAG vs Function heatmap (core modules)",
        "variance": "MAG vs Function heatmap (differential modules)",
    }.get(mode, "MAG vs Function heatmap")

    if representatives_only:
        title_text += " – representatives only"

    plt.suptitle(title_text, fontsize=13, y=1.06, fontweight="bold")

    # ---- Colorbars ----
    pos = ax.get_position()
    cax = fig.add_axes([pos.x0 - 0.06, pos.y0 - 0.23, 0.015, 0.14])
    cb = fig.colorbar(im, cax=cax, orientation="vertical")
    cb.ax.tick_params(labelsize=7)
    cb.set_label("Pathway completeness (%)", fontsize=8, labelpad=2)

    if size_im is not None:
        cax_size = fig.add_axes([pos.x1 + 0.04, pos.y0, 0.015, 0.14])
        cb_size = fig.colorbar(size_im, cax=cax_size, orientation="vertical")
        cb_size.ax.tick_params(labelsize=7)
        cb_size.set_label("MAG number", fontsize=8, labelpad=2)

    mod = "all" if top_modules is None else f"top{top_modules}"
    rep_tag = "representatives" if representatives_only else "all_mags"

    out_file = os.path.join(output_path, f"pathway_module_heatmap_{mode}_{rep_tag}_{mod}.{fmt}")
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[INFO] pathway module heatmap saved to {out_file}")
    return out_file