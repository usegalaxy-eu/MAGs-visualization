import pandas as pd
import numpy as np
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import time
import os
from .version import __version__
from .sanky_taxa import generate_taxa_sanky,taxa_sanky_rank
from .comp_conta_plot import completeness_contamination_plot, rank_completeness_contamination_plot, metadata_completeness_contamination_plot
from .heatmap import mag_heatmap
from .drep_cluster_plot import drep_cluster_plot
from .drep_cluster_func import drep_cluster_functional_plot
from .pathway_module_heatmap import pathway_module_heatmap
from .overview import OverviewOptions, PANELS, overview_figure


def positive_int(value):
    ivalue = int(value)
    if ivalue < 5:
        raise argparse.ArgumentTypeError(f"{value} must be >= 5")
    return ivalue


def int_or_none(value):
        if value.lower() == "none":
            return None
        return int(value)


def reset_matplotlib():
    """
    Reset matplotlib global state so plots do not inherit style/layout
    from previously created figures.
    (Required, if you want to use subcommand "all".)
    """
    plt.close("all")
    mpl.rcdefaults()
    plt.rcdefaults()


def run_clean(plot_func, *args, **kwargs):
    """
    Run a plotting function in a clean matplotlib state.
    """
    reset_matplotlib()
    return plot_func(*args, **kwargs)


CM_TO_INCH = 1 / 2.54

def cm(value):
    """Convert cm to inches for matplotlib."""
    return value * CM_TO_INCH

def _fig_size_tuple(args):
    """Convert user-supplied fig_size (in cm) to inches tuple."""
    if getattr(args, "fig_size", None) is None:
        return None
    width, height = args.fig_size
    return (cm(width), cm(height))


def read_table(path, index_col=None, prefer_tsv=None):
    if path is None:
        return None

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        first = fh.readline()

    # If caller explicitly prefers TSV
    if prefer_tsv is True:
        df = pd.read_csv(path, sep="\t", engine="python")
    elif prefer_tsv is False:
        # Prefer CSV
        try:
            df = pd.read_csv(path, sep=",", engine="python")
        except Exception:
            df = pd.read_csv(path, sep=None, engine="python")
    else:
        # Auto mode: tabs -> TSV else CSV
        if "\t" in first:
            df = pd.read_csv(path, sep="\t", engine="python")
        else:
            try:
                df = pd.read_csv(path, sep=",", engine="python")
            except Exception:
                df = pd.read_csv(path, sep=None, engine="python")

    # Only set index if requested
    if index_col is not None and df.shape[1] > index_col:
        df = df.set_index(df.columns[index_col])

    return df


def merged_coverm(coverm_dfs):
    """
    Merge multiple coverm tables column wise.
    """
    if not coverm_dfs:
        print("[WARN] No CoverM DataFrames to merge.")
        return pd.DataFrame()

    if len(coverm_dfs) == 1:
        df = coverm_dfs.get('coverm_0', next(iter(coverm_dfs.values())))
        return df

    dfs = [coverm_dfs[k] for k in sorted(coverm_dfs.keys())]
    coverm_merged = pd.concat(dfs, axis=1)
    print(f"[INFO] coverm merged: {coverm_merged.shape} rows x columns")
    return coverm_merged


def check_path(output_path):
    if output_path is None:
        raise SystemExit("[ERROR] --output is required.")
    if not os.path.exists(output_path):
        os.makedirs(output_path)
        print(f"[INFO] Created folder: {output_path}")
    else:
        print(f"[INFO] Folder already exists: {output_path}")


def resolve_col_names(metadata_df, col_indices):
    """
    Data_column returns 1-based index.
    """
    if col_indices is None:
        return None
    flat = []
    for c in col_indices:
        for part in str(c).split(","):
            part = part.strip()
            if part:
                flat.append(part)
    cols = metadata_df.columns.tolist()
    result = []
    for c in flat:
        try:
            idx = int(c) - 1
            result.append(cols[idx])
        except (ValueError, IndexError):
            result.append(c)
    return result

def is_index_input(meta_cols):
    if meta_cols is None:
        return False
    flat = []
    for c in meta_cols:
        for part in str(c).split(","):
            part = part.strip()
            if part:
                flat.append(part)
    return bool(flat) and all(p.isdigit() for p in flat)


# ---- Data loading ---- #
def load_dfs(coverm, checkm, checkm2, gtdb, drep, bakta=None, quast=None, metadata=None, pathways=None):
    """
    Load only what is provided.
    """
    dfs = {}

    # checkm
    if checkm is not None:
        dfs["checkm"] = read_table(checkm, index_col=0)
        print(f"[INFO] checkm loaded:  {dfs['checkm'].shape}")
    else:
        dfs["checkm"] = None

    # checkm2
    if checkm2 is not None:
        dfs["checkm2"] = read_table(checkm2, index_col=0)
        print(f"[INFO] checkm2 loaded: {dfs['checkm2'].shape}")
    else:
        dfs["checkm2"] = None

    # drep
    if drep is not None:
        dfs["drep"] = read_table(drep, index_col=None, prefer_tsv=None)
        print(f"[INFO] drep loaded:   {dfs['drep'].shape}")

        # dRep cleanup to ensure genome column exists
        dfs["drep"].columns = [str(c).strip() for c in dfs["drep"].columns]
        if "genome" not in dfs["drep"].columns and "Genome" in dfs["drep"].columns:
            dfs["drep"] = dfs["drep"].rename(columns={"Genome": "genome"})
        if "genome" not in dfs["drep"].columns:
            dfs["drep"] = dfs["drep"].reset_index().rename(columns={"index": "genome"})
    else:
        dfs["drep"] = None

    # gtdb
    if gtdb is not None:
        dfs["gtdb"] = read_table(gtdb, index_col=0, prefer_tsv=True)
        print(f"[INFO] gtdb loaded:   {dfs['gtdb'].shape}")
    else:
        dfs["gtdb"] = None

    # coverm (table or folder)
    coverm_dfs = {}
    if coverm is None:
        print("[INFO] No coverm provided.")
    else:
        if os.path.isdir(coverm):
            files = [f for f in os.listdir(coverm) if not f.startswith(".")]
            if not files:
                print(f"[WARN] No files found in coverm directory: {coverm}")
            for i, file in enumerate(sorted(files)):
                path = os.path.join(coverm, file)
                coverm_dfs[f"coverm_{i}"] = read_table(path, index_col=0)
                print(f"[INFO] coverm_{i} loaded: {coverm_dfs[f'coverm_{i}'].shape}")
        else:
            coverm_dfs["coverm_0"] = read_table(coverm, index_col=0)
            print(f"[INFO] coverm_0 loaded: {coverm_dfs['coverm_0'].shape}")

    dfs["coverm"] = coverm_dfs

    # optional tables
    dfs["bakta"] = read_table(bakta, index_col=0) if bakta is not None else None
    if dfs["bakta"] is not None:
        print(f"[INFO] bakta loaded:  {dfs['bakta'].shape}")

    dfs["quast"] = read_table(quast, index_col=0) if quast is not None else None
    if dfs["quast"] is not None:
        print(f"[INFO] quast loaded:  {dfs['quast'].shape}")

    dfs["metadata"] = read_table(metadata, index_col=None) if metadata is not None else None
    if dfs["metadata"] is not None:
        print(f"[INFO] metadata loaded:{dfs['metadata'].shape}")
    else:
        print("[INFO] No metadata file provided.")
    
    dfs["pathways"] = read_table(pathways, index_col=None, prefer_tsv=True) if pathways is not None else None
    if dfs["pathways"] is not None:
        print(f"[INFO] pathways loaded:{dfs['pathways'].shape}")
    return dfs


# Help texts of options shared by several subcommands. They are filled in by
# _add_default_help() wherever an option has no explicit help, so the meaning of
# an option is written once and identical in every subcommand and in docs/cli.md.
OPTION_HELP = {
    "--coverm": "CoverM relative abundance table (genomes x samples, incl. 'unmapped'), or a folder of tables",
    "--checkm": "CheckM quality table (Bin Id, Completeness, Contamination, ...)",
    "--checkm2": "CheckM2 quality report (Name, Completeness, Contamination, ...)",
    "--gtdb": "GTDB-Tk summary table (user_genome, classification, ...)",
    "--drep": "dRep cluster table Cdb.csv (genome, secondary_cluster, ...)",
    "--quast": "QUAST report with one column per MAG",
    "--bakta": "Bakta summary with one column per MAG",
    "--pathways": "KEGG pathway completeness table (contig, module_accession, completeness, pathway_name, "
                  "pathway_class)",
    "--metadata": "Sample metadata table; the first column holds the sample names used in the CoverM columns",
    "--meta_cols": "Metadata columns shown as annotation bars (names or 1-based column indices)",
    "--meta_col": "Metadata column used for colouring (name or 1-based column index)",
    "--meta_bin_width": "Bin width used to group numeric metadata values",
    "--tax_level": "GTDB rank used to group and colour the MAGs",
    "--tax_levels": "GTDB ranks shown as taxonomy columns next to the clusters",
    "--rank": "Rank of the second Sankey diagram (sankey_plot_rank_filtered), which links each MAG to its taxon "
              "at this rank",
    "--top_n": "Number of largest dRep clusters (or most frequent categories) shown; minimum 5",
    "--top_modules": "Number of KEGG modules shown ('none' for all)",
    "--output": "Output folder (created if it does not exist)",
    "--format": "Output file format",
    "--max_col": "Show the N most abundant taxa as columns and pool the rest as 'Others'",
    "--no_log": "Plot MAGs per taxon on a linear instead of a log10 axis",
    "--top_bar_height": "Height of the MAGs-per-taxon bar chart",
    "--hspace": "Vertical gap between the bar chart and the heatmap",
    "--heatmap_width": "Width of the heatmap",
    "--spacer_legend": "Gap between the legend and the metadata bars",
    "--spacer_meta": "Gap between the metadata bars and the heatmap",
    "--spacer_heatmap": "Gap between the heatmap and the MAGs-per-sample histogram",
    "--legend": "Width of the legend column",
    "--meta_bar_add": "Extra width added to the metadata bars",
    "--top_bar_spacer": "Gap above the MAGs-per-taxon bar chart",
    "--row_fontsize": "Font size (pt) of the MAG names on the rows",
    "--font_size": "Font size (pt) of node labels",
    "--quality": "Colour the completeness/contamination plot by MAG quality",
    "--tax": "Colour the completeness/contamination plot by taxonomy",
    "--color_by": "Colouring of the completeness/contamination plot",
    "--tax_levels_space": "Spacing between taxonomy columns (cm)",
    "--mode": "Module selection: mean = core modules, variance = differential modules",
    "--top_representatives": "Limit the heatmap to the representatives of the N largest clusters",
    "--column_choice": "Unused; accepted for backwards compatibility",
    "--show_module_labels": "Unused; accepted for backwards compatibility (see --module_label)",
}


def _add_default_help(parser):
    """Fill missing help texts from OPTION_HELP and show defaults in --help."""
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            for sp in action.choices.values():
                sp.formatter_class = argparse.ArgumentDefaultsHelpFormatter
                _add_default_help(sp)
        elif not action.help and action.option_strings:
            action.help = OPTION_HELP.get(action.option_strings[-1], None)


# ---- Parser builder - subcommands ---- #
def build_parser():
    parser = argparse.ArgumentParser(
        prog="mags-visualization",
        description="Visualize and analyse MAGs (Metagenome-Assembled Genomes).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-v", "-V", "--version", action="version", version=__version__)

    sub = parser.add_subparsers(dest="command", required=True)

    tax_levels = ["domain", "phylum", "class", "order", "family", "genus", "species"]

    # ---- sample-heatmap ----
    p = sub.add_parser("sample-heatmap", help="Create the sample heatmap plot.")
    p.add_argument("--coverm", required=True, dest="coverm_path")
    p.add_argument("--gtdb", required=True, dest="gtdb_file")
    p.add_argument("--metadata", dest="metadata_file", default=None)
    p.add_argument("--meta_cols", nargs="+", default=None)
    p.add_argument("--tax_level", choices=tax_levels, default="phylum")
    p.add_argument("-o", "--output", required=True, dest="output")
    p.add_argument("--format", choices=["png", "pdf", "svg"], default="png", dest="format")
    p.add_argument("--meta_bin_width", type=float, default=5.0, dest="meta_bin_width")
    p.add_argument("--max_col", type=int, default=10, dest="max_col")
    p.add_argument("--no_log", action="store_true")


    # layout options
    layout = p.add_argument_group("Heatmap Layout Options (all values in cm)")
    layout.add_argument("--top_bar_height", type=float, default=1.8, dest="top_bar_height")
    layout.add_argument("--hspace", type=float, default=0.6, dest="hspace")
    layout.add_argument("--heatmap_width", type=float, default=28.0, dest="heatmap_width")
    layout.add_argument("--spacer_legend", type=float, default=0.8, dest="spacer_legend")
    layout.add_argument("--spacer_meta", type=float, default=5.0, dest="spacer_meta")
    layout.add_argument("--spacer_heatmap", type=float, default=0.25, dest="spacer_heatmap")
    layout.add_argument("--legend", type=float, default=6.4, dest="legend")
    layout.add_argument("--meta_bar_add", type=float, default=3.8, dest="meta_bar_add")
    layout.add_argument("--top_bar_spacer", type=float, default=0.0, dest="top_bar_spacer")

    # ---- drep-cluster-annotation ----
    p = sub.add_parser("drep-cluster-annot", help="Create the dRep cluster plot.")
    p.add_argument("--drep", required=True, dest="drep_file")
    p.add_argument("--gtdb", required=True, dest="gtdb_file")
    p.add_argument("--checkm2", dest="checkm2_file", default=None)
    p.add_argument("--quast", dest="quast_file", default=None)
    p.add_argument("--bakta", dest="bakta_file", default=None)
    p.add_argument("--tax_levels", nargs="+", choices=tax_levels, default=["phylum", "genus"])
    p.add_argument("--top_n", type=positive_int, default=30, dest="n")
    p.add_argument("--tax_levels_space", type=float, default=0.8, dest="tax_levels_space", help="Spacing between taxonomy columns (cm)")
    p.add_argument("--fig_size", nargs=2, type=float, metavar=("WIDTH", "HEIGHT"), dest="fig_size",
                   help="Figure size in cm (width height)")
    p.add_argument("--require_quality", action="store_true", dest="require_quality",
                   help="Only show clusters whose representative has CheckM2 data")
    p.add_argument("-o", "--output", required=True, dest="output")
    p.add_argument("--format", choices=["png", "pdf", "svg"], default="png", dest="format")

    # ---- drep-cluster-functional ----
    p = sub.add_parser("drep-cluster-func", help="Create the dRep cluster plot with functional annotation.")
    p.add_argument("--drep", required=True, dest="drep_file")
    p.add_argument("--gtdb", required=True, dest="gtdb_file")
    p.add_argument("--checkm2", dest="checkm2_file", default=None)
    p.add_argument("--quast", dest="quast_file", default=None)
    p.add_argument("--bakta", dest="bakta_file", default=None)
    p.add_argument("--pathways", dest="pathways_file", default=None, help="TSV with contig, module accession, completeness, pathway_name, pathway_class, ...")
    p.add_argument("--tax_levels", nargs="+", choices=tax_levels, default=["phylum", "genus"])
    p.add_argument("--top_n", type=positive_int, default=30, dest="n")
    p.add_argument("--tax_levels_space", type=float, default=0.8, dest="tax_levels_space",
                   help="Spacing between taxonomy columns (cm)")
    p.add_argument("--fig_size", nargs=2, type=float, metavar=("WIDTH", "HEIGHT"), dest="fig_size",
                   help="Figure size in cm (width height)")
    p.add_argument("-o", "--output", required=True, dest="output")
    p.add_argument("--format", choices=["png" , "pdf", "svg"], default="png", dest="format")
    p.add_argument("--top_modules", type=int_or_none, default=35, dest="top_modules")
    p.add_argument("--mode", choices=["mean", "variance", "both"], default="mean", dest="mode",
                   help="mean = core modules, variance = differential modules, or both.")
    p.add_argument("--module_label", choices=["id", "name", "both"], default="id",
                   dest="module_label", help="Module labels on the x-axis: id, name or both")
    
    # ---- pathway-module-heatmap ----
    p = sub.add_parser("pathway-module-heatmap", help="Create a MAG vs pathway-module heatmap grouped by pathway_class.")
    p.add_argument("--pathways", required=True, dest="pathways_file",
                   help="TSV with contig, module_accession, completeness, pathway_name, pathway_class, ...")
    p.add_argument("--top_modules", type=int_or_none, default=30, dest="top_modules",
                   help="Number of modules to show. Use None to show all modules.")
    p.add_argument("--mode", choices=["mean", "variance"], default="mean", dest="mode",
                   help="mean = core modules, variance = differential modules")
    p.add_argument("--fig_size", nargs=2, type=float, metavar=("WIDTH", "HEIGHT"), dest="fig_size",
                   help="Figure size in cm (width height)")
    p.add_argument("--row_fontsize", type=int, default=8, dest="row_fontsize")
    p.add_argument("-o", "--output", required=True, dest="output")
    p.add_argument("--format", choices=["png", "pdf", "svg"], default="png", dest="format")
    p.add_argument("--drep", dest="drep_file", default=None, 
                   help="Optional dRep table used to restrict the heatmap to representatives only.")
    p.add_argument("--gtdb", dest="gtdb_file", default=None,
        help="Optional GTDB table required together with --drep when using --top_representatives.")
    p.add_argument("--top_representatives", type=int, default=None, dest="top_representatives",
               help="Limit number of representative MAGs (e.g. top 30)")
    p.add_argument("--sort_by", choices=["cluster_size", "none"], default="cluster_size", help="How to select top representatives")
    p.add_argument("--module_label", choices=["id", "name", "both"], default="id",
                   dest="module_label", help="Module labels on the x-axis: id, name, or both.")

    # ---- comp-conta ----
    p = sub.add_parser("comp-conta", help="Create completeness/contamination plots.")
    p.add_argument("--tools", choices=["both", "checkm", "checkm2"], default="both",
               help="Which tool results to plot (default: both)")
    p.add_argument("--checkm", required=False, dest="checkm_file")
    p.add_argument("--checkm2", required=False, dest="checkm2_file")
    p.add_argument("--gtdb", dest="gtdb_file", default=None)
    p.add_argument("--metadata", dest="metadata_file", default=None)
    p.add_argument("--meta_col", dest="meta_col", default=None)
    p.add_argument("--meta_bin_width", type=float, default=5.0, dest="meta_bin_width")
    p.add_argument("--top_n", type=positive_int, default=30, dest="n")
    p.add_argument("--mode", choices=["quality", "tax", "meta"], default="quality",
                   help="quality = plain plots; tax = rank plot needs --gtdb; meta = metadata plot needs --metadata and --meta_col")
    p.add_argument("--tax_level", choices=tax_levels, default="phylum", dest="tax_level")
    p.add_argument("--fig_size", nargs=2, type=float, metavar=("WIDTH", "HEIGHT"), dest="fig_size",
                   help="Figure size in cm (width height)")
    p.add_argument("-o", "--output", required=True, dest="output")
    p.add_argument("--format", choices=["png", "pdf", "svg"], default="png", dest="format")

    # ---- taxa-sankey ----
    p = sub.add_parser("taxa-sankey", help="Create GTDB taxa sankey plots.")
    p.add_argument("--gtdb", required=True, dest="gtdb_file")
    p.add_argument("--coverm", dest="coverm_path", default=None,
                   help="Optional CoverM abundance table; bar sizes become abundance-based")
    p.add_argument("--rank", choices=tax_levels, default="phylum", dest="rank")
    p.add_argument("--font_size", type=int, default=10, dest="font_size")
    p.add_argument("--order_links", action="store_true", dest="order_links",
                   help="Order nodes/links to reduce crossings (writes sankey_plot_ordered.*)")
    p.add_argument("--format", choices=["html", "png", "pdf", "svg"], default="html")
    p.add_argument("-o", "--output", required=True, dest="output")

    # ---- all (optional) ----
    # ---- overview-figure ----
    p = sub.add_parser(
        "overview-figure",
        help="Page-width multi-panel overview figure (abundance, quality, taxonomy, function).",
        description=(
            "Draw a publication-ready overview of a MAG catalogue at its final print size. "
            "Panels: abundance (relative abundance per sample and taxon, with sample metadata), "
            "quality (completeness vs contamination), taxonomy (flow diagram between two ranks) and "
            "function (largest dRep clusters with KEGG module completeness). Panels whose inputs are "
            "missing are skipped; use --panels to draw a subset or a single panel."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    g = p.add_argument_group("inputs")
    g.add_argument("--gtdb", required=True, dest="gtdb_file", help="GTDB-Tk summary table (all panels)")
    g.add_argument("--coverm", dest="coverm_path", default=None,
                   help="CoverM relative abundance table, or a folder of tables (abundance panel)")
    g.add_argument("--checkm2", dest="checkm2_file", default=None, help="CheckM2 quality report (quality panel)")
    g.add_argument("--checkm", dest="checkm_file", default=None,
                   help="CheckM quality report, used for the quality panel if no CheckM2 report is given")
    g.add_argument("--drep", dest="drep_file", default=None, help="dRep Cdb table (function panel)")
    g.add_argument("--pathways", dest="pathways_file", default=None,
                   help="KEGG pathway completeness table (function panel)")
    g.add_argument("--metadata", dest="metadata_file", default=None,
                   help="Sample metadata; the first column must hold the CoverM sample names (abundance panel)")
    g = p.add_argument_group("layout")
    g.add_argument("--panels", nargs="+", choices=PANELS, default=PANELS, help="Panels to draw, in figure order")
    g.add_argument("--width_cm", type=float, default=17.0, help="Figure width in cm (17 = full page width)")
    g.add_argument("--font_size", type=float, default=8.0,
                   help="Base font size in pt at the final figure size; dense labels use 85%% of it")
    g.add_argument("--color_rank", choices=tax_levels, default="phylum",
                   help="Rank whose taxa get one colour each, shared by all panels")
    g = p.add_argument_group("abundance panel")
    g.add_argument("--tax_level", choices=tax_levels, default="phylum", help="Rank of the heatmap rows")
    g.add_argument("--max_taxa", type=int, default=10, help="Show the N most abundant taxa, pool the rest as Other")
    g.add_argument("--meta_cols", nargs="+", default=None,
                   help="Metadata columns shown as rows above the heatmap (names or 1-based indices). "
                        "Boolean columns become filled/empty squares, others coloured squares")
    g.add_argument("--meta_labels", nargs="+", default=None, help="Display names for --meta_cols (same order)")
    g.add_argument("--sample_label_col", default=None,
                   help="Metadata column used as sample label (name or 1-based index)")
    g.add_argument("--sample_group_col", default=None,
                   help="Metadata column grouping the samples (sorted together, separated by a gap); "
                        "name or 1-based index")
    g.add_argument("--group_order", nargs="+", default=None, help="Order of the sample groups")
    g = p.add_argument_group("quality panel")
    g.add_argument("--hq_completeness", type=float, default=90.0, help="Minimum completeness of high-quality MAGs")
    g.add_argument("--hq_contamination", type=float, default=5.0, help="Maximum contamination of high-quality MAGs")
    g.add_argument("--mq_completeness", type=float, default=70.0, help="Minimum completeness of medium-quality MAGs")
    g.add_argument("--mq_contamination", type=float, default=None,
                   help="Maximum contamination of medium-quality MAGs (no limit if unset)")
    g.add_argument("--max_contamination_axis", type=float, default=8.0,
                   help="Upper limit of the contamination axis; MAGs above it are drawn as triangles at the edge "
                        "and their values listed in the legend")
    g = p.add_argument_group("taxonomy panel")
    g.add_argument("--flow_ranks", nargs=2, choices=tax_levels, default=["phylum", "genus"],
                   metavar=("FROM", "TO"), help="First and last rank of the flow diagram")
    g = p.add_argument_group("function panel")
    g.add_argument("--top_n", type=int, default=30, help="Number of largest dRep clusters shown")
    g.add_argument("--top_modules", type=int, default=35, help="Number of KEGG modules shown")
    g.add_argument("--module_mode", choices=["mean", "variance"], default="mean",
                   help="Select modules with the highest mean completeness (core functions) or variance "
                        "(differences between clusters)")
    g = p.add_argument_group("output")
    g.add_argument("--format", nargs="+", choices=["png", "pdf", "svg"], default=["png"],
                   help="One or more output formats")
    g.add_argument("--dpi", type=int, default=300, help="Resolution of PNG output")
    g.add_argument("--prefix", default="overview_figure", help="File name (without extension) inside --output")
    g.add_argument("-o", "--output", required=True, dest="output", help="Output folder")

    p = sub.add_parser("all", help="Run all plots (legacy behaviour).")
    p.add_argument("--coverm", dest="coverm_path", default=None)
    p.add_argument("--checkm", dest="checkm_file", default=None)
    p.add_argument("--checkm2", dest="checkm2_file", default=None)
    p.add_argument("--gtdb", dest="gtdb_file", default=None)
    p.add_argument("--drep", dest="drep_file", default=None)
    p.add_argument("--quast", dest="quast_file", default=None)
    p.add_argument("--bakta", dest="bakta_file", default=None)
    p.add_argument("--metadata", dest="metadata_file", default=None)
    p.add_argument("--pathways", dest="pathways_file", default=None)
    p.add_argument("--rank", choices=tax_levels, default="phylum", dest="rank")
    p.add_argument("--tax_level", choices=tax_levels, default=None, dest="tax_level")
    p.add_argument("--tax_levels", nargs="+", choices=tax_levels, default=["phylum", "genus"])
    p.add_argument("--top_n", type=positive_int, default=30, dest="n")
    p.add_argument("--quality", action="store_true")
    p.add_argument("--require_quality", action="store_true", dest="require_quality",
                   help="Only show clusters whose representative has CheckM2 data")
    p.add_argument("--tax", action="store_true")
    p.add_argument("--meta_col", dest="meta_col", default=None)
    p.add_argument("--meta_cols", nargs="+", default=None)
    p.add_argument("--meta_bin_width", type=float, default=5.0, dest="meta_bin_width")
    p.add_argument("--column_choice", nargs="+", metavar="METRIC", dest="column_choice")
    p.add_argument("--color_by", choices=["quality", "tax", "meta"], dest="color_by")
    p.add_argument("--tax_levels_space", type=float, default=0.8, dest="tax_levels_space")
    p.add_argument("--fig_size", nargs=2, type=float, metavar=("WIDTH", "HEIGHT"), dest="fig_size",
                   help="Figure size in cm (width height)")
    p.add_argument("-o", "--output", required=True, dest="output")
    p.add_argument("--format", choices=["png", "pdf", "svg"], default="png", dest="format")
    # heatmap layout subset
    p.add_argument("--top_bar_height", type=float, default=1.8, dest="top_bar_height")
    p.add_argument("--hspace", type=float, default=0.6, dest="hspace")
    p.add_argument("--heatmap_width", type=float, default=28.0, dest="heatmap_width")
    p.add_argument("--spacer_legend", type=float, default=0.8, dest="spacer_legend")
    p.add_argument("--spacer_meta", type=float, default=5.0, dest="spacer_meta")
    p.add_argument("--spacer_heatmap", type=float, default=0.25, dest="spacer_heatmap")
    p.add_argument("--legend", type=float, default=6.4, dest="legend")
    p.add_argument("--meta_bar_add", type=float, default=3.8, dest="meta_bar_add")
    p.add_argument("--top_bar_spacer", type=float, default=0.0, dest="top_bar_spacer")
    p.add_argument("--max_col", type=int, default=10, dest="max_col")
    p.add_argument("--no_log", action="store_true")
    # functional annoation
    p.add_argument("--top_modules", type=int_or_none, default=35, dest="top_modules")
    p.add_argument("--mode", choices=["mean", "variance"], default="mean", dest="mode")
    p.add_argument("--show_module_labels", action="store_true", dest="show_module_labels")
    p.add_argument("--row_fontsize", type=int, default=8, dest="row_fontsize")
    p.add_argument("--top_representatives", type=int, default=None, dest="top_representatives")

    _add_default_help(parser)
    return parser


def parse_arguments(argv=None):
    return build_parser().parse_args(argv)


# ---- Subcommand implementation ---- #
def run_sample_heatmap(args):
    dfs = load_dfs(
        coverm=args.coverm_path,
        checkm=None,
        checkm2=None,
        gtdb=args.gtdb_file,
        drep=None,
        metadata=args.metadata_file,
    )
    dfs["coverm"] = merged_coverm(dfs["coverm"])
    check_path(args.output)

    # Mode: name or index
    meta_cols = args.meta_cols
    if dfs.get("metadata") is not None and meta_cols is not None:
        if is_index_input(meta_cols):
            meta_cols = resolve_col_names(dfs["metadata"], meta_cols)
            print(f"[INFO] Resolved meta_cols by index: {meta_cols}")
        else:
            print(f"[INFO] Using meta_cols by name: {meta_cols}")

    if dfs.get("metadata") is not None:
        first_col = dfs["metadata"].columns[0]
        dfs["metadata"] = dfs["metadata"].set_index(first_col)
        print(f"[INFO] Metadata index set to: '{first_col}'")

    mag_heatmap(
        dfs["coverm"],
        dfs["gtdb"],
        args.output,
        rank=args.tax_level,
        metadata_df=dfs.get("metadata"),
        meta_cols=meta_cols,
        meta_bin_width=args.meta_bin_width,
        fmt=args.format,
        top_bar_height=cm(args.top_bar_height),
        hspace=cm(args.hspace),
        heatmap_width=cm(args.heatmap_width),
        spacer_legend=cm(args.spacer_legend),
        spacer_meta=cm(args.spacer_meta),
        spacer_heatmap=cm(args.spacer_heatmap),
        legend=cm(args.legend),
        meta_bar_add=cm(args.meta_bar_add),
        top_bar_spacer=cm(args.top_bar_spacer),
        max_col=args.max_col,
        log_top=not args.no_log,
    )


def run_drep_cluster(args):
    dfs = load_dfs(
        coverm=None,
        checkm=None,
        checkm2=args.checkm2_file,
        gtdb=args.gtdb_file,
        drep=args.drep_file,
        bakta=args.bakta_file,
        quast=args.quast_file,
        metadata=None,
    )
    check_path(args.output)

    drep_cluster_plot(
        dfs["drep"],
        dfs["gtdb"],
        args.output,
        tax_levels=args.tax_levels,
        top_n=args.n,
        fig_size=_fig_size_tuple(args),
        fmt=args.format,
        tax_levels_space=cm(args.tax_levels_space),
        checkm2_df=dfs.get("checkm2"),
        quast_df=dfs.get("quast"),
        bakta_df=dfs.get("bakta"),
        require_quality=getattr(args, "require_quality", False),
    )


def run_drep_cluster_func(args):
    dfs = load_dfs(
        coverm=None,
        checkm=None,
        checkm2=args.checkm2_file,
        gtdb=args.gtdb_file,
        drep=args.drep_file,
        pathways=args.pathways_file,
        bakta=args.bakta_file,
        quast=args.quast_file,
        metadata=None,
    )
    check_path(args.output)

    drep_cluster_functional_plot(
        dfs["drep"],
        dfs["gtdb"],
        dfs["pathways"],
        args.output,
        tax_levels=args.tax_levels,
        top_n=args.n,
        top_modules=args.top_modules,
        fig_size=_fig_size_tuple(args),
        fmt=args.format,
        mode=args.mode,
        tax_levels_space=cm(args.tax_levels_space),
        module_label=args.module_label,
    )


def run_pathway_module_heatmap(args):
    dfs = load_dfs(
        coverm=None,
        checkm=None,
        checkm2=None,
        gtdb=args.gtdb_file,
        drep=args.drep_file,
        bakta=None,
        quast=None,
        metadata=None,
        pathways=args.pathways_file,
    )
    check_path(args.output)

    pathway_module_heatmap(
        pathway_df=dfs["pathways"],
        output_path=args.output,
        top_modules=args.top_modules,
        mode=args.mode,
        fmt=args.format,
        fig_size=_fig_size_tuple(args),
        row_fontsize=args.row_fontsize,
        module_label=args.module_label,
        representatives_df=dfs.get("drep"),
        gtdb_df=dfs.get("gtdb"),
        top_representatives=args.top_representatives,
    )


def run_comp_conta(args):
    # Validate that required files are provided for the selected tools
    if args.tools in ("checkm", "both") and args.checkm_file is None:
        raise SystemExit("[ERROR] --checkm is required when --tools is 'checkm' or 'both'")
    if args.tools in ("checkm2", "both") and args.checkm2_file is None:
        raise SystemExit("[ERROR] --checkm2 is required when --tools is 'checkm2' or 'both'")

    dfs = load_dfs(
        coverm=None,
        checkm=args.checkm_file,
        checkm2=args.checkm2_file,
        gtdb=args.gtdb_file,
        drep=None,
        metadata=args.metadata_file,
    )
    check_path(args.output)

    # Mode: name or index
    meta_col = args.meta_col
    if dfs.get("metadata") is not None and meta_col is not None:
        if is_index_input([meta_col]):
            resolved = resolve_col_names(dfs["metadata"], [meta_col])
            meta_col = resolved[0] if resolved else meta_col
            print(f"[INFO] Resolved meta_col by index: {meta_col}")
        else:
            print(f"[INFO] Using meta_cols by name: {meta_col}")
        
    if dfs.get("metadata") is not None:
        first_col = dfs["metadata"].columns[0]
        dfs["metadata"] = dfs["metadata"].set_index(first_col)


    fig_size = _fig_size_tuple(args) or (9, 8)

    if args.mode == "quality":
        if args.tools in ("checkm", "both"):
            completeness_contamination_plot(dfs["checkm"], args.output, tag="checkm",
                                            title="CheckM: Completeness vs Contamination",
                                            fig_size=fig_size, fmt=args.format)
        if args.tools in ("checkm2", "both"):
            completeness_contamination_plot(dfs["checkm2"], args.output, tag="checkm2",
                                            title="CheckM2: Completeness vs Contamination",
                                            fig_size=fig_size, fmt=args.format)
        return

    if args.mode == "tax":
        if dfs["gtdb"] is None:
            raise SystemExit("[ERROR] comp-conta --mode tax requires --gtdb")
        rank_completeness_contamination_plot(
            dfs["checkm"],
            dfs["checkm2"],
            dfs["gtdb"],
            args.tax_level,
            args.output,
            args.n,
            fig_size=_fig_size_tuple(args) or (10, 8),
            fmt=args.format,
            tools=args.tools,
        )
        return

    if args.mode == "meta":
        if dfs["metadata"] is None or meta_col is None:
            raise SystemExit("[ERROR] comp-conta --mode meta requires --metadata and --meta_col")
        if meta_col not in dfs["metadata"].columns:
            raise SystemExit(f"[ERROR] meta_col '{meta_col}' not found. Available: {list(dfs['metadata'].columns)}")
        metadata_completeness_contamination_plot(
            dfs["checkm"],
            dfs["checkm2"],
            dfs["metadata"],
            meta_col=meta_col,
            output_path=args.output,
            fig_size=_fig_size_tuple(args) or (10, 8),
            fmt=args.format,
            bin_width=args.meta_bin_width,
            top_n=args.n,
            tools=args.tools,
        )
        return

    raise SystemExit(f"[ERROR] Unknown comp-conta mode: {args.mode}")


def run_taxa_sankey(args):
    dfs = load_dfs(
        coverm=args.coverm_path,
        checkm=None,
        checkm2=None,
        gtdb=args.gtdb_file,
        drep=None,
        metadata=None,
    )
    check_path(args.output)
    coverm_df = merged_coverm(dfs["coverm"]) if dfs["coverm"] else None
    generate_taxa_sanky(dfs["gtdb"], args.output, args.rank, fmt=args.format,
                        coverm=coverm_df, font_size=args.font_size,
                        order_links=getattr(args, "order_links", False))
    taxa_sanky_rank(dfs["gtdb"], args.output, args.rank, fmt=args.format)


def run_overview_figure(args):
    dfs = load_dfs(
        coverm=args.coverm_path,
        checkm=args.checkm_file if args.checkm2_file is None else None,
        checkm2=args.checkm2_file,
        gtdb=args.gtdb_file,
        drep=args.drep_file,
        metadata=args.metadata_file,
        pathways=args.pathways_file,
    )
    check_path(args.output)
    meta_cols = args.meta_cols
    label_col, group_col = args.sample_label_col, args.sample_group_col
    if dfs.get("metadata") is not None:
        # Galaxy data_column parameters pass 1-based indices
        if is_index_input(meta_cols):
            meta_cols = resolve_col_names(dfs["metadata"], meta_cols)
        if label_col is not None and is_index_input([label_col]):
            label_col = resolve_col_names(dfs["metadata"], [label_col])[0]
        if group_col is not None and is_index_input([group_col]):
            group_col = resolve_col_names(dfs["metadata"], [group_col])[0]
    try:
        opt = OverviewOptions(
            panels=args.panels, width_cm=args.width_cm, font_size=args.font_size,
            abundance_rank=args.tax_level, max_taxa=args.max_taxa, meta_cols=meta_cols,
            meta_labels=args.meta_labels, sample_label_col=label_col,
            sample_group_col=group_col, group_order=args.group_order,
            hq_completeness=args.hq_completeness, hq_contamination=args.hq_contamination,
            mq_completeness=args.mq_completeness, mq_contamination=args.mq_contamination,
            max_contamination_axis=args.max_contamination_axis, flow_ranks=tuple(args.flow_ranks),
            top_n=args.top_n, top_modules=args.top_modules, module_mode=args.module_mode,
            color_rank=args.color_rank,
        )
        coverm = merged_coverm(dfs["coverm"]) if dfs.get("coverm") else None
        quality = dfs.get("checkm2") if dfs.get("checkm2") is not None else dfs.get("checkm")
        reset_matplotlib()
        overview_figure(
            os.path.join(args.output, args.prefix), opt, coverm=coverm, gtdb=dfs["gtdb"], quality=quality,
            drep=dfs.get("drep"), pathways=dfs.get("pathways"), metadata=dfs.get("metadata"),
            formats=args.format, dpi=args.dpi,
        )
    except ValueError as e:
        raise SystemExit(f"[ERROR] {e}")


def run_all(args):
    """
    Runs all based on provided inputs.
    """
    dfs = load_dfs(
        args.coverm_path,
        args.checkm_file,
        args.checkm2_file,
        args.gtdb_file,
        args.drep_file,
        bakta=args.bakta_file,
        quast=args.quast_file,
        metadata=args.metadata_file,
        pathways=args.pathways_file
    )
    check_path(args.output)

    if dfs.get("coverm"):
        dfs["coverm"] = merged_coverm(dfs["coverm"])

    comp_fig_size = _fig_size_tuple(args)

    tax_rank = args.tax_level or args.rank or "phylum"

    # sankey
    if dfs["gtdb"] is not None:
        reset_matplotlib()
        generate_taxa_sanky(dfs["gtdb"], args.output, args.rank)

        reset_matplotlib()
        taxa_sanky_rank(dfs["gtdb"], args.output, args.rank)

    # comp/conta
    run_quality = args.quality or (not args.quality and not args.tax)
    run_tax = args.tax or (not args.quality and not args.tax)

    if run_quality and dfs["checkm"] is not None:
        run_clean(
            completeness_contamination_plot,
            dfs["checkm"], args.output, tag="checkm",
            title="CheckM: Completeness vs Contamination",
            fig_size=comp_fig_size or (9, 8), fmt=args.format
        )
    if run_quality and dfs["checkm2"] is not None:
        run_clean(
            completeness_contamination_plot,
            dfs["checkm2"], args.output, tag="checkm2",
            title="CheckM2: Completeness vs Contamination",
            fig_size=comp_fig_size or (9, 8), fmt=args.format
        )

    if run_tax and dfs["gtdb"] is not None and dfs["checkm"] is not None and dfs["checkm2"] is not None:
        run_clean(
            rank_completeness_contamination_plot,
            dfs["checkm"], dfs["checkm2"], dfs["gtdb"], tax_rank,
            args.output, args.n,
            fig_size=comp_fig_size or (10, 8), fmt=args.format
        )

    # drep-cluster-annot
    if dfs.get("drep") is not None and dfs.get("gtdb") is not None:
        run_clean(
            drep_cluster_plot,
            dfs["drep"], dfs["gtdb"], args.output,
            tax_levels=args.tax_levels, top_n=args.n,
            fig_size=comp_fig_size, fmt=args.format,
            tax_levels_space=cm(args.tax_levels_space),
            checkm2_df=dfs.get("checkm2"),
            quast_df=dfs.get("quast"),
            bakta_df=dfs.get("bakta"),
            require_quality=getattr(args, "require_quality", False),
       )
    
    # drep-cluster-func
    if dfs.get("drep") is not None and dfs.get("gtdb") is not None and dfs.get("pathways") is not None:
        run_clean(
            drep_cluster_functional_plot,
            dfs["drep"],
            dfs["gtdb"],
            dfs["pathways"],
            args.output,
            top_n=args.n,
            fmt=args.format,
            tax_levels=args.tax_levels,
            top_modules=args.top_modules,
            fig_size=_fig_size_tuple(args),
            tax_levels_space=cm(args.tax_levels_space),
        )
    
    # pathway-heatmap-completeness
    if dfs.get("pathways") is not None:
        run_clean(
            pathway_module_heatmap,
            pathway_df=dfs["pathways"],
            output_path=args.output,
            top_modules=args.top_modules,
            mode=args.mode,
            fmt=args.format,
            fig_size=_fig_size_tuple(args),
            row_fontsize=args.row_fontsize,
            representatives_df=dfs.get("drep"),
            gtdb_df=dfs.get("gtdb"),
            top_representatives=args.top_representatives,
        )

    # heatmap
    if dfs.get("coverm") is not None and dfs.get("gtdb") is not None:
        run_clean(
            mag_heatmap,
            dfs["coverm"], dfs["gtdb"], args.output,
            rank=tax_rank,
            metadata_df=dfs.get("metadata"),
            meta_cols=args.meta_cols,
            meta_bin_width=args.meta_bin_width,
            fmt=args.format,
            top_bar_height=cm(args.top_bar_height),
            hspace=cm(args.hspace),
            heatmap_width=cm(args.heatmap_width),
            spacer_legend=cm(args.spacer_legend),
            spacer_meta=cm(args.spacer_meta),
            spacer_heatmap=cm(args.spacer_heatmap),
            legend=cm(args.legend),
            meta_bar_add=cm(args.meta_bar_add),
            top_bar_spacer=cm(args.top_bar_spacer),
            max_col=args.max_col,
            log_top=not args.no_log,
        )


# ---- Main ---- #
def main(argv=None):
    start_time = time.time()
    args = parse_arguments(argv)

    if args.command == "sample-heatmap":
        run_sample_heatmap(args)
    elif args.command == "drep-cluster-annot":
        run_drep_cluster(args)
    elif args.command == "drep-cluster-func":
        run_drep_cluster_func(args)
    elif args.command == "pathway-module-heatmap":
        run_pathway_module_heatmap(args)
    elif args.command == "comp-conta":
        run_comp_conta(args)
    elif args.command == "taxa-sankey":
        run_taxa_sankey(args)
    elif args.command == "overview-figure":
        run_overview_figure(args)
    elif args.command == "all":
        run_all(args)
    else:
        raise SystemExit(f"[ERROR] Unknown command: {args.command}")

    end_time = time.time()
    print(f"[INFO] Run time: {time.strftime('%H:%M:%S', time.gmtime(end_time - start_time))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
