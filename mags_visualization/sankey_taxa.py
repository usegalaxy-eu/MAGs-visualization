import pandas as pd
import numpy as np
import plotly.graph_objects as go
import os
import re

from mags_visualization.id_normalizer import normalize_genome_id

prefix_map = {
    "d": "domain",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species"
}

def normalize_rank(rank: str) -> str:
    """Accepts short codes (d/p/c/o/f/g/s) or full names; returns full rank name."""
    if not isinstance(rank, str):
        raise ValueError("rank must be a string")
    r = rank.strip().lower()
    # map short code -> full
    short2full = {
        "d": "domain", "p": "phylum", "c": "class",
        "o": "order",  "f": "family", "g": "genus", "s": "species"
    }
    return short2full.get(r, r)  # if already 'phylum', returns 'phylum'


def clean_tax_label(label: str, rank_name: str) -> str:
    """
    Convert gtdb labels ("g_Gilliamella") into display labels.
    Missing values = unkown.
    """
    if pd.isna(label) or not str(label).strip():
        return f"Unknown {rank_name.capitalize()}"
    
    label = str(label).strip()
    cleaned = re.sub(r"^[a-z]__*", "", label).strip()
    
    if cleaned:
        return cleaned
    return f"Unknown {rank_name.capitalize()}"


def is_missing_taxon(value) -> bool:
    """
    Return True if a taxonomy entry is missing or unclassified
    """
    if pd.isna(value):
        return True

    value = str(value).strip()
    if not value:
        return True

    # Remove GTDB prefix like d__, p__, g__
    cleaned = re.sub(r"^[a-z]__*", "", value).strip()

    if cleaned == "":
        return True

    cleaned_lower = cleaned.lower()

    # Treat generic unknown/unclassified labels as missing
    if cleaned_lower.startswith("unclassified"):
        return True
    if cleaned_lower.startswith("unknown"):
        return True
    if cleaned_lower in {"uncultured", "unidentified", "na", "n/a", "none"}:
        return True

    return False

def _rendered_link_y_keys(node_order, links_df):
    """
    Render a minimal sankey (same node values, pad, thickness and size as the
    final figure) and return, for each link in array order, the rendered
    vertical position (source mid y, target mid y) extracted from the SVG.
    Returns None if the number of rendered links does not match the input.
    """
    node_map = {n: i for i, n in enumerate(node_order)}
    tmp = links_df.copy()
    tmp["source_idx"] = tmp["source"].map(node_map)
    tmp["target_idx"] = tmp["target"].map(node_map)

    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=25,
            thickness=20,
            line=dict(color="black", width=0),
            label=[str(n) for n in node_order],
            color=["rgba(0,0,0,0)"] * len(node_order)
        ),
        link=dict(
            source=tmp["source_idx"],
            target=tmp["target_idx"],
            value=tmp["value"]
        )
    )])
    fig.update_layout(width=2000, height=1600)

    svg = fig.to_image(format="svg").decode("utf-8")
    paths = re.findall(r'<path class="sankey-link"[^>]*d="([^"]+)"', svg)
    if len(paths) != len(links_df):
        return None

    keys = []
    for d in paths:
        nums = [float(v) for v in re.findall(r"[-.\de+]+", d)]
        sy_top, sy_bot = nums[1], nums[15]
        ty_top, ty_bot = nums[7], nums[9]
        keys.append(((sy_top + sy_bot) / 2, (ty_top + ty_bot) / 2))
    return keys


def _minimize_crossings(links_df):
    """
    Reorder nodes within each column and re-sort links so that d3-sankey
    produces as few crossings as possible.

    d3-sankey stacks the ribbons leaving a node strictly in array order, so
    the only way to keep them parallel is to sort the links by their final
    rendered vertical position. Node positions are independent of the link
    order, so we bootstrap: render a temporary sankey, read the rendered
    (source y, target y) of every link back out of the SVG, sort the links by
    those coordinates, and repeat until the order is stable.

    Nodes are ordered per column with the barycenter heuristic (median
    position of neighbours in the adjacent column) first.
    """
    layers = ["domain", "phylum", "class", "order", "family", "genus", "species", "_total"]

    all_nodes = pd.Index(pd.concat([links_df["source"], links_df["target"]]).unique())

    def layer_of(node):
        if node == "_total":
            return len(layers) - 1
        return layers.index(str(node).split("|")[0])

    pred = {n: set() for n in all_nodes}
    succ = {n: set() for n in all_nodes}
    for _, row in links_df.iterrows():
        succ[row["source"]].add(row["target"])
        pred[row["target"]].add(row["source"])

    columns = [[] for _ in layers]
    for n in all_nodes:
        columns[layer_of(n)].append(n)

    # barycenter sweeps (median of neighbour positions in the adjacent column)
    for _ in range(10):
        # left -> right, order each column by its predecessors
        for L in range(1, len(layers)):
            pos_prev = {n: i for i, n in enumerate(columns[L - 1])}

            def key(n):
                ns = [pos_prev[m] for m in pred[n] if m in pos_prev]
                if ns:
                    return (float(np.median(ns)), 0)
                ns2 = [pos_prev[m] for m in succ[n] if m in pos_prev]
                return (float(np.median(ns2)) if ns2 else 0.0, 1)

            columns[L] = sorted(columns[L], key=key)

        # right -> left, order each column by its successors
        for L in range(len(layers) - 2, -1, -1):
            pos_next = {n: i for i, n in enumerate(columns[L + 1])}

            def key(n):
                ns = [pos_next[m] for m in succ[n] if m in pos_next]
                if ns:
                    return (float(np.median(ns)), 0)
                ns2 = [pos_next[m] for m in pred[n] if m in pos_next]
                return (float(np.median(ns2)) if ns2 else 0.0, 1)

            columns[L] = sorted(columns[L], key=key)

    ordered_nodes = [n for col in columns for n in col]
    node_map = {name: i for i, name in enumerate(ordered_nodes)}

    links_df = links_df.copy()
    links_df["source_idx"] = links_df["source"].map(node_map)
    links_df["target_idx"] = links_df["target"].map(node_map)

    # bootstrap link sort on the rendered vertical positions
    for _ in range(4):
        keys = _rendered_link_y_keys(ordered_nodes, links_df)
        if keys is None:
            break
        order = sorted(range(len(keys)), key=lambda i: keys[i])
        reordered = links_df.iloc[order]
        if reordered["source_idx"].tolist() == links_df["source_idx"].tolist() and \
           reordered["target_idx"].tolist() == links_df["target_idx"].tolist():
            links_df = reordered
            break
        links_df = reordered.reset_index(drop=True)

    return ordered_nodes, links_df


def generate_taxa_sanky(gtdb, output_path, rank, fmt="html", coverm=None, font_size=13, order_links=False):
    if rank is not None:
        rank = normalize_rank(rank)

    tax_columns = ["domain", "phylum", "class", "order", "family", "genus", "species"]

    # optional per-MAG abundance from a CoverM table (summed over samples);
    # when given, link values (and thus bar sizes) are abundance-based
    mag_abundance = {}
    if coverm is not None:
        cov = coverm.copy()
        if cov.index.name == "Genome":
            cov = cov.reset_index()
        elif "Genome" not in cov.columns:
            cov = cov.reset_index()
            cov = cov.rename(columns={cov.columns[0]: "Genome"})
        sample_cols = cov.select_dtypes(include="number").columns
        cov["user_genome"] = cov["Genome"].map(normalize_genome_id)
        mag_abundance = dict(zip(cov["user_genome"], cov[sample_cols].sum(axis=1)))

    tax_split = gtdb.reset_index()["classification"].str.split(";", expand=True)
    tax_split.columns = tax_columns[:tax_split.shape[1]]
    tax_split = tax_split.reindex(columns=tax_columns)
    tax_split = tax_split.replace({"": None, " ": None})
    tax_split["user_genome"] = gtdb.reset_index().iloc[:, 0]

    if rank is not None and rank not in tax_columns:
        raise ValueError(f"Requested rank '{rank}' is not in parsed columns {tax_columns}")

    links = []

    for _, row in tax_split.iterrows():
        if mag_abundance:
            value = mag_abundance.get(normalize_genome_id(row["user_genome"]), 0.0)
        else:
            value = 1.0

        valid_ranks = []

        # collect only real classified ranks
        for current_rank in tax_columns:
            current_value = row[current_rank]
            if not is_missing_taxon(current_value):
                valid_ranks.append((current_rank, str(current_value).strip()))

        # skip completely unclassified rows
        if not valid_ranks:
            continue

        # build links between known ranks
        previous_node = None
        for current_rank, current_value in valid_ranks:
            current_node = f"{current_rank}|{current_value}"

            if previous_node is not None:
                links.append((previous_node, current_node, value))

            previous_node = current_node

        # if species is missing, add one terminal Unknown Species node
        species_value = row["species"]
        if is_missing_taxon(species_value):
            last_rank, last_value = valid_ranks[-1]

            if last_rank == "genus":
                unknown_species_node = f"species|unknown_from|{last_rank}|{last_value}"
                links.append((previous_node, unknown_species_node, value))

    if not links:
        raise ValueError("No classified taxonomy paths were found for Sankey plotting.")

    links_df = (
        pd.DataFrame(links, columns=["source", "target", "value"])
        .groupby(["source", "target"], as_index=False)["value"].sum()
    )

    # invisible "total" sink column right of species: plotly draws labels on
    # the LEFT for the rightmost column, so species needs a real column after
    # it (with links, since unlinked nodes are dropped) to keep its labels
    # on the right
    total_node = "_total"
    sink_links = []
    species_nodes = links_df.loc[
        links_df["target"].str.startswith("species|"), "target"
    ].unique()
    for species_node in species_nodes:
        value = links_df.loc[links_df["target"] == species_node, "value"].sum()
        sink_links.append({"source": species_node, "target": total_node, "value": value})

    if sink_links:
        links_df = pd.concat([
            links_df,
            pd.DataFrame(sink_links)
        ], ignore_index=True)

    if order_links:
        nodes, links_df = _minimize_crossings(links_df)
    else:
        nodes = pd.Index(pd.concat([links_df["source"], links_df["target"]]).unique())
        node_map = {name: i for i, name in enumerate(nodes)}
        links_df = links_df.copy()
        links_df["source_idx"] = links_df["source"].map(node_map)
        links_df["target_idx"] = links_df["target"].map(node_map)

    rank_colors = {
        "domain": "#1f77b4",
        "phylum": "#ff7f0e",
        "class": "#2ca02c",
        "order": "#d62728",
        "family": "#9467bd",
        "genus": "#8c564b",
        "species": "#e377c2"
    }

    clean_labels = []
    node_colors = []

    for node in nodes:
        if node == total_node:
            node_colors.append("rgba(0,0,0,0)")
            clean_labels.append("")
            continue

        parts = str(node).split("|")
        node_rank = parts[0]
        node_colors.append(rank_colors.get(node_rank, "lightgray"))

        # normal node: "rank|value"
        if len(parts) == 2:
            raw_value = parts[1]
            clean_labels.append(clean_tax_label(raw_value, node_rank))

        # unknown species: "species|unknown_from|genus|g__Gilliamella"
        elif len(parts) == 4 and parts[1] == "unknown_from":
            parent_rank = parts[2]
            parent_value = parts[3]
            parent_label = clean_tax_label(parent_value, parent_rank)
            clean_labels.append(f"Unknown Species")

        else:
            clean_labels.append(node)

    link_colors = [
        "rgba(0,0,0,0)" if target == total_node else "rgba(0,0,0,0.3)"
        for target in links_df["target"]
    ]

    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=25,
            thickness=20,
            # per-node line arrays are ignored by plotly's sankey trace, so a
            # transparent sink node would still get a visible outline; width=0
            # removes the outline globally (the sink fill stays transparent)
            line=dict(color="black", width=0),
            label=clean_labels,
            color=node_colors
        ),
        link=dict(
            source=links_df["source_idx"],
            target=links_df["target_idx"],
            value=links_df["value"],
            color=link_colors
        )
    )])

    fig.update_layout(
        title_text="Taxonomic Classification Sankey",
        font=dict(size=font_size),
        width=2000,
        height=1600
    )

    os.makedirs(output_path, exist_ok=True)
    stem = "sankey_plot_ordered" if order_links else "sankey_plot"
    out_file = os.path.join(output_path, f"{stem}.{fmt}")

    if fmt == "html":
        fig.write_html(out_file)
    else:
        fig.write_image(out_file)


def taxa_sanky_rank(gtdb, output_path, rank, fmt="html"):

    tax_split = gtdb.reset_index()["classification"].str.split(";", expand=True)
    tax_split.columns = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    tax_split = tax_split.replace({"": None, " ": None})

    df_sankey = gtdb.reset_index()[["user_genome"]].copy()
    df_sankey[rank] = tax_split[rank]

    links_df = df_sankey.groupby(["user_genome", rank]).size().reset_index(name="count")
    links_df.columns = ["source", "target", "count"]

    nodes = pd.Index(pd.concat([links_df["source"], links_df["target"]]).unique())
    node_map = {name: i for i, name in enumerate(nodes)}
    links_df["source_idx"] = links_df["source"].map(node_map)
    links_df["target_idx"] = links_df["target"].map(node_map)

    rank_colors = {
        "genome": "#7f7f7f",  
        "domain": "#1f77b4",
        "phylum": "#ff7f0e",
        "class": "#2ca02c",
        "order": "#d62728",
        "family": "#9467bd",
        "genus": "#8c564b",
        "species": "#e377c2"
    }

    node_colors = []
    for node in nodes:
        if node in df_sankey["user_genome"].values:
            node_colors.append(rank_colors["genome"])
        else:
            node_colors.append(rank_colors.get(rank, "lightgray"))

    clean_labels = []
    for label in nodes:
        label_str = str(label)

        if label_str.endswith(".fasta"):
            clean_labels.append(label_str)
        elif "__" in label_str:
            val = label_str.split("__")[-1]
            clean_labels.append(val if val else f"Unknown {rank.capitalize()}")
        else:
            clean_labels.append(label_str)

    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=15,
            line=dict(color="black", width=0.5),
            label=clean_labels,
            color=node_colors
        ),
        link=dict(
            source=links_df["source_idx"],
            target=links_df["target_idx"],
            value=links_df["count"]
        )
    )])

    fig.update_layout(
        title_text=f"Genome → {rank.capitalize()} Sankey",
        font=dict(size=10),
        width=1200,
        height=800
    )

    out_file = os.path.join(output_path, f"sankey_plot_rank_filtered.{fmt}")

    if fmt == "html":
        fig.write_html(out_file)
    else:
        fig.write_image(out_file)