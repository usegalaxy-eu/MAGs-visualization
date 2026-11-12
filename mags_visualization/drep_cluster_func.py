import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
from matplotlib.colors import Normalize
from .id_normalizer import normalize_genome_id, format_genome_display
from .drep_cluster_plot import extract_rank


def simplify_pathway_class(value: str) -> str:
    """
    Return second semicolon separated name of kegg pathway class.
    """
    if pd.isna(value):
        return "Unclassified"

    parts = [p.strip() for p in str(value).split(";") if str(p).strip()]
    if len(parts) >= 2:
        return parts[1]
    if len(parts) == 1:
        return parts[0]
    return "Unclassified"


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


def prepare_cluster_data(drep_df: pd.DataFrame, gtdb_df: pd.DataFrame, tax_levels, top_n: int, pathway_df=None):
    """Build table of top_n dRep cluster with their representatives and taxonomy."""
    drep = drep_df.copy()
    gtdb = gtdb_df.copy()

    if drep.index.name is not None:
        drep = drep.reset_index()

    if "genome" not in drep or "secondary_cluster" not in drep:
        raise ValueError("dRep must contain 'genome' and 'secondary_cluster'")
    
    drep["genome_normalized"] = drep["genome"].apply(normalize_genome_id)

    # Count genomes and keep top n biggest clusters
    cluster_counts = (drep.groupby("secondary_cluster").size().reset_index(name="n_members")
        .sort_values("n_members", ascending=False).head(top_n))

    total_clusters = drep["secondary_cluster"].nunique()

    if gtdb.index.name is not None:
        gtdb = gtdb.reset_index()

    if "user_genome" not in gtdb.columns or "classification" not in gtdb.columns:
        raise ValueError("GTDB file must contain 'user_genome' and 'classification' columns")

    gtdb["user_genome_normalized"] = gtdb["user_genome"].apply(normalize_genome_id)

    for lvl in set(tax_levels):
        gtdb[lvl] = gtdb["classification"].apply(lambda x: extract_rank(x, lvl))

    gtdb_genomes = set(gtdb["user_genome_normalized"].values)

    pw_genomes = set()
    if pathway_df is not None:
        pw_genomes = set(pathway_df["contig"].apply(normalize_genome_id).unique())

    print(f"[INFO] Creating dRep cluster plot (top {top_n}, tax_levels={tax_levels})")

    def get_representative(cluster_genomes):
        # Prefer genomes present in pathway and gtdb
        if pw_genomes:
            for genome in cluster_genomes:
                if genome in gtdb_genomes and genome in pw_genomes:
                    return genome
        for genome in cluster_genomes:
            if genome in gtdb_genomes:
                return genome
        return cluster_genomes.iloc[0]

    cluster_reps = (drep.groupby("secondary_cluster")["genome_normalized"].apply(get_representative).reset_index())
    cluster_reps.columns = ["secondary_cluster", "representative"]

    # Combine cluster counts with representative genomes
    cluster_data = cluster_counts.merge(cluster_reps, on="secondary_cluster", how="left")
    cluster_data["representative_display"] = (cluster_data["representative"].apply(format_genome_display))

    cols_to_merge = ["user_genome_normalized"] + list(tax_levels)
    cluster_data = cluster_data.merge(
        gtdb[cols_to_merge],
        left_on="representative",
        right_on="user_genome_normalized",
        how="left",)

    for lvl in tax_levels:
        cluster_data[lvl] = cluster_data.get(lvl, "Unclassified").fillna("Unclassified")

    # sort for plotting
    cluster_data = cluster_data.sort_values("n_members").reset_index(drop=True)

    return cluster_data, drep, total_clusters

# ---- Module Selection ---- #
def select_top_modules(pw: pd.DataFrame, representatives, top_modules: int | None = 35, mode="mean"):
    """
    Filter pathway data ranked by mean completeness or variance across representatives.
    
    Select modules based on:
    - 'mean' = core functions
    - 'variance' = differences between clusters
    """
    if top_modules is None:
        return pw

    # Keep only pathway rows for selected representative genomes
    pw = pw[pw["contig_normalized"].isin(representatives)]

    # Build matrix
    # rows = genomes, columns = modules, values = completeness
    pivot = pw.pivot_table(
        index="contig_normalized",
        columns="module_accession",
        values="completeness",
        aggfunc="max"
    )
    
    pivot = pivot.loc[representatives]

    # Score each module depending on selected mode
    if mode == "mean":
        scores = pivot.fillna(0).mean(axis=0)
    elif mode == "variance":
        scores = pivot.fillna(0).var(axis=0)
    else:
        raise ValueError("mode must be 'mean' or 'variance'")

    top = scores.sort_values(ascending=False).head(top_modules).index 
    return pw[pw["module_accession"].isin(top)].copy()


def prepare_pathway_matrix(pathway_df: pd.DataFrame, representatives: list[str], top_modules: int | None = 35, mode="mean",):
    """Create heatmap matrix: 
    rows = representative MAGs
    columns = modules accessions
    values = completeness"""

    pw = pathway_df.copy()

    required_cols = {"contig", "module_accession", "completeness", "pathway_name", "pathway_class",}
    missing = required_cols - set(pw.columns)
    if missing:
        raise ValueError(f"Pathway file is missing required columns: {sorted(missing)}")

    pw["contig_normalized"] = pw["contig"].apply(normalize_genome_id)
    pw["completeness"] = pd.to_numeric(pw["completeness"], errors="coerce")
    pw["pathway_class_simple"] = pw["pathway_class"].apply(simplify_pathway_class)

    pw = pw[pw["contig_normalized"].isin(representatives)].copy()

    if pw.empty:
        raise ValueError("No pathway rows matched the selected cluster representatives")

    pw = select_top_modules(pw, representatives=representatives, top_modules=top_modules, mode=mode)

    module_meta = (
        pw[["module_accession", "pathway_name", "pathway_class_simple"]]
        .drop_duplicates()
        .sort_values(["pathway_class_simple", "module_accession"])
        .reset_index(drop=True))

    module_order = module_meta["module_accession"].tolist()

    heatmap_df = (pw.pivot_table(index="contig_normalized", columns="module_accession",
            values="completeness", aggfunc="max").reindex(index=representatives, columns=module_order))
    
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

    return heatmap_df, module_meta, class_blocks


def plot_single_functional(cluster_data: pd.DataFrame, drep: pd.DataFrame, total_clusters: int,
    heatmap_df: pd.DataFrame, module_meta: pd.DataFrame, class_blocks, output_path: str,
    tax_levels, top_n: int, fmt: str, suffix: str, tax_levels_space: float = 0.3, 
    fig_size=None, module_label: str = "id"):
    """
    Create a single functional cluster plot:
    - horizontal bar chart (MAG counts per cluster)
    - taxonomy annotation
    - kegg pathway completeness heatmap
    """

    if heatmap_df.empty:
        raise ValueError(f"Heatmap matrix is empty for plot '{suffix}'")

    if fig_size is None:
        heatmap_width = max(8, len(module_meta) * 0.22)
        fig_width = 6 + len(tax_levels) * 1.6 + heatmap_width
        fig_height = max(8, top_n * 0.32)
        fig_size = (fig_width, fig_height)

    fig = plt.figure(figsize=fig_size)
    n_tax = len(tax_levels)

    # [Sample labels | bar plot | spacer | taxonomy | spacer | heatmap]
    width_ratios = ([0.08, 4.8, 0.02] + [max(1.1, tax_levels_space * 3.5)] * n_tax
        + [0.35, max(5.5, len(module_meta) * 0.20)])

    gs = gridspec.GridSpec(
        1,
        len(width_ratios),
        figure=fig,
        width_ratios=width_ratios,
        wspace=0.04,
    )

    ax_samples = fig.add_subplot(gs[0, 0])
    ax_bars = fig.add_subplot(gs[0, 1])

    tax_axes = []
    tax_start_col = 3
    for i, lvl in enumerate(tax_levels):
        tax_axes.append((lvl, fig.add_subplot(gs[0, tax_start_col + i])))

    ax_heatmap = fig.add_subplot(gs[0, len(width_ratios) - 1])

    y_pos = np.arange(len(cluster_data))

    # ---- Bar plot: number of MAGs per cluster ----
    ax_bars.barh(
        y_pos,
        cluster_data["n_members"].values,
        height=0.7,
        color="#1f2937",
        edgecolor="#111827",
        linewidth=0.5,
    )

    # X-axis scale and tick spacing
    max_members = int(cluster_data["n_members"].max())
    nice_max = int(np.ceil(max_members * 1.08))
    step = 5 if nice_max <= 50 else 10

    ax_bars.set_xlim(0, nice_max)
    ax_bars.set_xticks(np.arange(0, nice_max + 1, step))
    ax_bars.set_ylim(-0.5, len(cluster_data) - 0.5)
    ax_bars.set_xlabel("# of MAGs in species-level cluster", fontsize=11)
    ax_bars.grid(axis="x", color="#e5e7eb", linewidth=0.5, zorder=0)
    ax_bars.set_axisbelow(True)
    ax_bars.set_yticks([])
    ax_bars.spines["left"].set_visible(False)
    ax_bars.spines["top"].set_visible(False)
    ax_bars.spines["right"].set_visible(False)

    # Count labels
    for i, (_, row) in enumerate(cluster_data.iterrows()):
        count = int(row["n_members"])
        ax_bars.text(
            count + nice_max * 0.01,
            i,
            str(count),
            va="center",
            ha="left",
            fontsize=9,
            fontweight="bold",
        )

    # ---- Representative genomes ----
    ax_samples.set_xlim(0, 1)
    ax_samples.set_ylim(-0.5, len(cluster_data) - 0.5)
    ax_samples.set_yticks(y_pos)
    ax_samples.set_yticklabels(cluster_data["representative_display"].values, fontsize=8, ha="right")
    ax_samples.tick_params(axis="y", which="both", length=10, width=0.4, pad=4)
    ax_samples.set_xticks([])
    for spine in ax_samples.spines.values():
        spine.set_visible(False)

    # ---- Taxonomy annotation ----
    for lvl, ax in tax_axes:
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.5, len(cluster_data) - 0.5)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(cluster_data[lvl].values, fontsize=8, ha="left")
        ax.set_xticks([])
        ax.tick_params(axis="y", which="both", length=0, pad=2)

        for spine in ax.spines.values():
            spine.set_visible(False)

        # Taxonomic title
        ax.text(
            0.0,
            1.02,
            lvl.capitalize(),
            ha="left",
            va="bottom",
            fontsize=9,
            fontweight="bold",
            transform=ax.transAxes,
        )

    # ---- Pathway completeness heatmap ----
    heatmap_display_df = heatmap_df.iloc[::-1] # flipping, so largest cluster at the top

    values = heatmap_display_df.to_numpy(dtype=float)
    masked = np.ma.masked_invalid(values)

    cmap = plt.get_cmap("Reds").copy()
    cmap.set_bad(color="white")
    norm = Normalize(vmin=0, vmax=100)

    im = ax_heatmap.imshow(
        masked,
        aspect="auto",
        interpolation="none",
        cmap=cmap,
        norm=norm,
    )

    n_rows, n_cols = values.shape

    ax_heatmap.set_xlim(-0.5, n_cols - 0.5)
    ax_heatmap.set_ylim(n_rows - 0.5, -0.5)
    ax_heatmap.set_yticks([])

    ax_heatmap.set_xticks(np.arange(n_cols))

    module_labels = build_module_labels(module_meta, module_label=module_label)

    ax_heatmap.set_xticklabels(
        module_labels,
        rotation=90,
        fontsize=8 if module_label in {"name", "both"} or n_cols > 25 else 9,
    )
    ax_heatmap.tick_params(axis="x", bottom=True, labelbottom=True, top=False, length=2)

    # Add grid lines between cells
    ax_heatmap.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax_heatmap.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax_heatmap.grid(which="minor", color="#d1d5db", linewidth=0.35)
    ax_heatmap.tick_params(which="minor", bottom=False, left=False)

    for spine in ax_heatmap.spines.values():
        spine.set_visible(False)

    # Separator lines between pathway classes
    for _, _, end in class_blocks[:-1]:
        ax_heatmap.axvline(end + 0.5, color="white", linewidth=2)

    # Colors for pathway class grouping labels
    class_palette = [
        "#4C78A8", "#F58518", "#54A24B", "#E45756", "#72B7B2",
        "#B279A2", "#FF9DA6", "#9D755D", "#593122", "#2E91E5"
    ]

    # Colored bar directly above the heatmap
    y_bar = -0.95
    bar_height = 0.22

    # text above connector line
    label_levels = [-1.55, -2.25, -2.95, -3.65]
    line_bottom_pad = -0.01  # gap between line and label

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
        ax_heatmap.add_patch(rect)

        y_text = label_levels[i % len(label_levels)]

        ax_heatmap.plot(
            [center, center],
            [y_bar, y_text - line_bottom_pad],
            color=color,
            linewidth=1.0,
            clip_on=False,
            zorder=5,
        )

        # Class label text
        ax_heatmap.text(
            center,
            y_text,
            cls,
            ha="center",
            va="bottom",
            fontsize=8,
            rotation=0,
            color=color,
            clip_on=False,
            zorder=6,
        )

    # ---- Colorbar ----
    pos = ax_heatmap.get_position()
    cax = fig.add_axes([pos.x0 - 0.06, pos.y0 - 0.20, 0.02, 0.12])
    cb = fig.colorbar(im, cax=cax, orientation="vertical")
    cb.ax.tick_params(labelsize=7)
    cb.set_label("Pathway completeness (%)", fontsize=7, labelpad=2)

    # ---- Info box ----
    total_mags = len(drep)
    mags_in_top_clusters = drep[drep["secondary_cluster"].isin(cluster_data["secondary_cluster"])].shape[0]

    info_text = (
        f"Total MAGs: {total_mags}\n"
        f"Total clusters: {total_clusters}\n"
        f"MAGs in top {top_n}: {mags_in_top_clusters} "
        f"({mags_in_top_clusters / total_mags * 100:.1f}%)"
    )

    pos_samples = ax_samples.get_position()

    fig.text(
        pos_samples.x0,
        pos_samples.y0 - 0.08,
        info_text,
        fontsize=10,
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3, pad=0.5),
        family='monospace'
    )

    title_text = {"core": "core functional modules", "difference": "differential functional modules",}.get(suffix, suffix)

    plt.suptitle(f"Taxonomic and {title_text.lower()} of the {top_n} largest dRep clusters ",
        fontsize=13, y=1.06, fontweight="bold",)

    out_file = os.path.join(output_path,f"drep_cluster_functional_{suffix}_top{top_n}.{fmt}")
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[INFO] dRep cluster functional plot saved to {out_file}")
    return out_file


# ---- Main plot function ----
def drep_cluster_functional_plot(drep_df: pd.DataFrame, gtdb_df: pd.DataFrame, pathway_df: pd.DataFrame,
                     output_path: str, tax_levels=("phylum", "genus"), top_n: int = 30, top_modules: int | None = 35, 
                     fmt: str = "png", mode: str = "mean", tax_levels_space: float = 0.3, 
                     fig_size=None, module_label: str ="id"):
    """
    Create a dRep cluster plot with a functional kegg pathway heatmap.

    Left: top clusters by MAG count
    Middle: taxonomy columns
    Right: pathway completeness heatmap
           rows = representative MAGs
           columns = KEGG module_accession
           grouped by pathway_class
    """
    
    os.makedirs(output_path, exist_ok=True)

    # ---- Validate tax levels ----
    valid_ranks = {"domain", "phylum", "class", "order", "family", "genus", "species"}
    
    if tax_levels is None:
        tax_levels = ("phylum", "genus")

    if isinstance(tax_levels, str):
        tax_levels = tuple(tax_levels.split())

    tax_levels = tuple([t.lower().strip() for t in tax_levels if str(t).strip() != ""])
    if len(tax_levels) == 0:
        tax_levels = ("phylum", "genus")

    for lvl in tax_levels:
        if lvl not in valid_ranks:
            raise ValueError(f"Invalid tax level '{lvl}'. Choose from {sorted(valid_ranks)}.")

    # ---- Choose mode ----
    mode = str(mode).lower().strip()
    valid_modes = {"mean", "variance", "both"}
    if mode not in valid_modes:
        raise ValueError(f"False mode '{mode}'. Choose from {valid_modes}.")

    print(f"[INFO] Creating dRep cluster plot (top {top_n}, tax_levels={tax_levels})")

    # Collect genomes that have pathway data
    pw_genomes = set(pathway_df["contig"].apply(normalize_genome_id).unique())
    print(f"[INFO] {len(pw_genomes)} genomes found in pathway file")

    # Fetch more clusters than needed, because some will be filtered out later
    fetch_n = top_n * 2
    cluster_data, drep, total_clusters = prepare_cluster_data(
        drep_df=drep_df,
        gtdb_df=gtdb_df,
        tax_levels=tax_levels,
        top_n=fetch_n,
    )

    # Remove clusters with no pathway data
    before = len(cluster_data)
    cluster_data = cluster_data[
        cluster_data["representative"].isin(pw_genomes)].reset_index(drop=True)
    after = len(cluster_data)

    if before != after:
        print(f"[WARN] {before - after} clusters removed (representative has no pathway data). "
              f"{after} clusters remaining.")

    # cluster data sorted - tail() to keep largest clusters
    cluster_data = cluster_data.tail(top_n).reset_index(drop=True)

    if len(cluster_data) < top_n:
        print(f"[WARN] Only {len(cluster_data)} clusters available after filtering (requested {top_n})")

    representatives = cluster_data["representative"].tolist()

    outputs = []

    # map user mode
    plot_modes = []
    if mode in {"mean", "both"}:
        plot_modes.append(("mean", "core"))
    if mode in {"variance", "both"}:
        plot_modes.append(("variance", "difference"))
    
    for select_mode, suffix in plot_modes:
        heatmap_df_core, module_meta_core, class_blocks_core = prepare_pathway_matrix(
            pathway_df=pathway_df,
            representatives=representatives,
            top_modules=top_modules,
            mode=select_mode,
        )

        out_file = plot_single_functional(
            cluster_data=cluster_data,
            drep=drep,
            total_clusters=total_clusters,
            heatmap_df=heatmap_df_core,
            module_meta=module_meta_core,
            class_blocks=class_blocks_core,
            output_path=output_path,
            tax_levels=tax_levels,
            top_n=top_n,
            fmt=fmt,
            suffix=suffix,
            tax_levels_space=tax_levels_space,
            fig_size=fig_size,
            module_label=module_label,
        )
        outputs.append(out_file)

    return outputs
