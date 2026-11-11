"""Page-width overview figure of a MAG catalogue.

The figure combines up to four panels, each of which can also be drawn on its own:

* ``abundance`` (A): relative abundance per sample and taxon, with sample metadata,
  MAGs per taxon and MAGs detected per sample.
* ``quality`` (B): completeness against contamination with marginal histograms.
* ``taxonomy`` (C): taxonomic composition as a flow (Sankey) diagram between ranks.
* ``function`` (D): the largest species-level (dRep secondary) clusters with the KEGG
  module completeness of their representative MAG.

All sizes are given in centimetres and font sizes in points, so the figure can be
made at its final print size (e.g. 17 cm for a full-width journal figure) and the
text is legible as printed. The layout is measured after a first drawing pass so
tick labels of one panel never run into the next panel.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, PathPatch, Rectangle
from matplotlib.path import Path as MPath

from .drep_cluster_func import prepare_cluster_data, prepare_pathway_matrix
from .id_normalizer import normalize_genome_id

CM = 1 / 2.54
PT_CM = 2.54 / 72  # one typographic point in cm
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
PANELS = ["abundance", "quality", "taxonomy", "function"]
PANEL_LETTERS = dict(zip(PANELS, "ABCD"))

# Colour-blind safe categorical palette (Okabe-Ito, then Tableau for overflow)
CATEGORICAL = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9",
               "#F0E442", "#000000", "#8C564B", "#7F7F7F", "#BCBD22", "#17BECF"]
OTHER_COLOR = "#B0B0B0"
DARK = "#333333"
HQ_COLOR, MQ_COLOR, LQ_COLOR = "#B2182B", "#8C8C8C", "#67A9CF"
CLASS_COLORS = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860",
                "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD"]
CLASS_SHORT = {
    "Amino acid metabolism": "Amino acid",
    "Carbohydrate metabolism": "Carbohydrate",
    "Energy metabolism": "Energy",
    "Lipid metabolism": "Lipid",
    "Metabolism of cofactors and vitamins": "Cofactors &\nvitamins",
    "Nucleotide metabolism": "Nucleotide",
    "Biosynthesis of other secondary metabolites": "Secondary\nmetabolites",
    "Biosynthesis of terpenoids and polyketides": "Terpenoids &\npolyketides",
    "Glycan metabolism": "Glycan",
    "Metabolism of other amino acids": "Other amino\nacids",
    "Xenobiotics biodegradation": "Xenobiotics",
    "Gene set": "Gene set",
}
TRUE_VALUES = {"true", "yes", "y", "1", "t"}
FALSE_VALUES = {"false", "no", "n", "0", "f", ""}


# --------------------------------------------------------------------------- #
# Options
# --------------------------------------------------------------------------- #
@dataclass
class OverviewOptions:
    """Settings of the overview figure (CLI flags map 1:1 onto these fields)."""

    panels: list[str] = field(default_factory=lambda: list(PANELS))
    width_cm: float = 17.0
    font_size: float = 8.0
    # abundance panel
    abundance_rank: str = "phylum"
    max_taxa: int = 10
    meta_cols: list[str] | None = None
    meta_labels: list[str] | None = None
    sample_label_col: str | None = None
    sample_group_col: str | None = None
    group_order: list[str] | None = None
    # quality panel
    quality_tool: str = "checkm2"
    hq_completeness: float = 90.0
    hq_contamination: float = 5.0
    mq_completeness: float = 70.0
    mq_contamination: float | None = None
    max_contamination_axis: float = 8.0
    # taxonomy panel
    flow_ranks: tuple[str, str] = ("phylum", "genus")
    # function panel
    top_n: int = 30
    top_modules: int = 35
    module_mode: str = "mean"
    # colouring (rank whose taxa get the categorical colours in all panels)
    color_rank: str = "phylum"

    def __post_init__(self):
        unknown = [p for p in self.panels if p not in PANELS]
        if unknown:
            raise ValueError(f"Unknown panel(s) {unknown}; choose from {PANELS}")
        for r in (self.abundance_rank, self.color_rank, *self.flow_ranks):
            if r not in RANKS:
                raise ValueError(f"Unknown rank {r!r}; choose from {RANKS}")
        if self.meta_labels and len(self.meta_labels) != len(self.meta_cols or []):
            raise ValueError("meta_labels needs one label per meta_cols entry")
        if RANKS.index(self.flow_ranks[0]) >= RANKS.index(self.flow_ranks[1]):
            raise ValueError("flow_ranks must go from a higher to a lower rank, e.g. phylum genus")
        if self.module_mode not in ("mean", "variance"):
            raise ValueError("module_mode must be 'mean' or 'variance'")

    # derived text sizes (pt)
    @property
    def fs(self) -> float:
        return self.font_size

    @property
    def fs_small(self) -> float:
        return round(self.font_size * 0.85, 2)

    @property
    def fs_letter(self) -> float:
        return round(self.font_size * 1.5, 2)

    def text_h_cm(self, size: float | None = None) -> float:
        """Height of one line of text in cm."""
        return (size or self.fs) * PT_CM * 1.2


# --------------------------------------------------------------------------- #
# Data preparation (pure functions, unit tested)
# --------------------------------------------------------------------------- #
def taxonomy_table(gtdb: pd.DataFrame) -> pd.DataFrame:
    """GTDB-Tk summary -> one row per (normalised) genome, one column per rank.

    Missing ranks (``g__`` etc.) become ``"Unclassified"``.
    """
    df = gtdb.reset_index() if "user_genome" not in gtdb.columns else gtdb
    rows = {}
    for genome, classification in zip(df["user_genome"], df["classification"]):
        parts = str(classification).split(";") if pd.notna(classification) else []
        ranks = {}
        for rank, part in zip(RANKS, parts + [""] * len(RANKS)):
            name = part.strip()[3:] if len(part.strip()) > 3 else ""
            ranks[rank] = name or "Unclassified"
        rows[normalize_genome_id(genome)] = ranks
    return pd.DataFrame.from_dict(rows, orient="index", columns=RANKS)


def coverm_relative_abundance(coverm: pd.DataFrame) -> pd.DataFrame:
    """CoverM table -> genomes x samples (% of reads), the ``unmapped`` row dropped.

    Accepts the genome either as index or as a ``Genome`` column and strips the
    `` Relative Abundance (%)`` suffix CoverM adds to sample columns.
    """
    df = coverm.copy()
    if "Genome" in df.columns:
        df = df.set_index("Genome")
    df.index = df.index.astype(str)
    df = df[df.index.str.lower() != "unmapped"]
    df.columns = [str(c).split(" Relative Abundance")[0].split("_R1")[0].strip() for c in df.columns]
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    df.index = df.index.map(normalize_genome_id)
    return df


def abundance_by_taxon(abund: pd.DataFrame, tax: pd.DataFrame, rank: str, max_taxa: int | None = 10):
    """Sum MAG abundance per sample and taxon.

    Returns ``(heat, mags_per_taxon, mags_per_sample)`` where ``heat`` is
    taxa x samples (taxa ordered by total abundance, the rest pooled as
    ``Other``), ``mags_per_taxon`` counts MAGs per taxon and ``mags_per_sample``
    the MAGs with abundance > 0 per sample.
    """
    taxon = tax.reindex(abund.index)[rank].fillna("Unclassified")
    heat = abund.groupby(taxon).sum()
    counts = taxon.value_counts()
    order = heat.sum(axis=1).sort_values(ascending=False).index.tolist()
    order = [t for t in order if t != "Unclassified"] + (["Unclassified"] if "Unclassified" in order else [])
    if max_taxa and len([t for t in order if t != "Unclassified"]) > max_taxa:
        keep = [t for t in order if t != "Unclassified"][:max_taxa]
        rest = [t for t in order if t not in keep and t != "Unclassified"]
        heat.loc["Other"] = heat.loc[rest].sum()
        counts.loc["Other"] = counts.loc[rest].sum()
        order = keep + ["Other"] + (["Unclassified"] if "Unclassified" in order else [])
    heat = heat.loc[order]
    mags_per_taxon = counts.reindex(order).fillna(0).astype(int)
    mags_per_sample = (abund > 0).sum(axis=0).astype(int)
    return heat, mags_per_taxon, mags_per_sample


def quality_category(completeness, contamination, hq_completeness=90.0, hq_contamination=5.0,
                     mq_completeness=70.0, mq_contamination=None) -> np.ndarray:
    """Classify MAGs as ``HQ``, ``MQ`` or ``LQ`` (everything else)."""
    comp = np.asarray(completeness, dtype=float)
    cont = np.asarray(contamination, dtype=float)
    hq = (comp >= hq_completeness) & (cont <= hq_contamination)
    mq = comp >= mq_completeness
    if mq_contamination is not None:
        mq &= cont <= mq_contamination
    return np.where(hq, "HQ", np.where(mq, "MQ", "LQ"))


def taxonomy_flow(tax: pd.DataFrame, ranks: list[str]):
    """Tree of taxa between ``ranks`` for a flow diagram.

    Nodes are identified by their lineage tuple so the flow stays a tree even if
    a name occurs under two parents. Returns ``(columns, values, children)``:
    ``columns[k]`` lists the nodes of rank ``k`` in drawing order (children follow
    their parent, largest first, ``Unclassified`` last), ``values[node]`` the
    number of genomes and ``children[node]`` the ordered child nodes.
    """
    paths = tax[ranks].astype(str)
    values: dict[tuple, int] = {}
    for k in range(len(ranks)):
        for key, n in paths.groupby(ranks[: k + 1]).size().items():
            values[key if isinstance(key, tuple) else (key,)] = int(n)

    def sort_key(node):
        return (node[-1] == "Unclassified", -values[node], node[-1])

    children: dict[tuple, list[tuple]] = {}
    for node in values:
        if len(node) < len(ranks):
            kids = [n for n in values if len(n) == len(node) + 1 and n[:-1] == node]
            children[node] = sorted(kids, key=sort_key)
    columns: list[list[tuple]] = [[] for _ in ranks]

    def visit(node):
        columns[len(node) - 1].append(node)
        for child in children.get(node, []):
            visit(child)

    for root in sorted([n for n in values if len(n) == 1], key=sort_key):
        visit(root)
    return columns, values, children


def spread_labels(y: np.ndarray, min_gap: float, lo: float, hi: float) -> np.ndarray:
    """Move 1-D label positions apart by at least ``min_gap`` within [lo, hi].

    Labels keep their order; overlapping runs are centred on their mean target.
    """
    y = np.asarray(y, dtype=float)
    order = np.argsort(y, kind="stable")
    target = y[order]
    pos = target.copy()
    for _ in range(200):
        moved = False
        # build clusters of labels that are too close and centre them
        start = 0
        while start < len(pos):
            end = start
            while end + 1 < len(pos) and pos[end + 1] - pos[end] < min_gap - 1e-9:
                end += 1
            if end > start:
                centre = target[start:end + 1].mean()
                span = (end - start) * min_gap
                new = centre - span / 2 + np.arange(end - start + 1) * min_gap
                if not np.allclose(new, pos[start:end + 1]):
                    pos[start:end + 1] = new
                    moved = True
            start = end + 1
        # keep inside bounds
        if pos[0] < lo:
            pos += lo - pos[0]
        if pos[-1] > hi:
            pos -= pos[-1] - hi
        if not moved:
            break
    out = np.empty_like(pos)
    out[order] = pos
    return out


def taxon_colors(taxa, palette=CATEGORICAL) -> dict[str, str]:
    """Stable colour per taxon; ``Other``/``Unclassified`` are grey."""
    colors, i = {}, 0
    for t in taxa:
        if t in ("Other", "Unclassified") or t in colors:
            colors.setdefault(t, OTHER_COLOR)
            continue
        colors[t] = palette[i % len(palette)]
        i += 1
    return colors


def to_bool(series: pd.Series) -> pd.Series | None:
    """Return a boolean series if every value looks boolean, else ``None``."""
    s = series.astype(str).str.strip().str.lower()
    if s.isin(TRUE_VALUES | FALSE_VALUES).all():
        return s.isin(TRUE_VALUES)
    return None


def sample_metadata(metadata: pd.DataFrame | None, samples: list[str], label_col=None, group_col=None,
                    group_order=None) -> pd.DataFrame:
    """Align metadata (first column = sample id) to the CoverM samples.

    Adds ``_label`` (display label), ``_group`` and returns rows in plotting
    order (by group, in ``group_order`` if given, then by label).
    """
    if metadata is None:
        meta = pd.DataFrame(index=pd.Index(samples, name="sample"))
    else:
        meta = metadata.copy()
        meta = meta.set_index(meta.columns[0]) if meta.index.name is None or meta.index.dtype != object else meta
        meta.index = meta.index.astype(str)
        meta = meta.reindex(samples)
    for col in (label_col, group_col):
        if col is not None and col not in meta.columns:
            raise ValueError(f"Metadata column {col!r} not found; available: {list(meta.columns)}")
    meta["_label"] = meta[label_col].astype(str) if label_col else meta.index
    meta["_label"] = meta["_label"].where(meta["_label"].notna() & (meta["_label"] != "nan"), meta.index.to_series())
    meta["_group"] = meta[group_col].astype(str) if group_col else ""
    if group_order:
        rank = {g: i for i, g in enumerate(group_order)}
        meta["_group_rank"] = meta["_group"].map(rank).fillna(len(rank))
    else:
        meta["_group_rank"] = pd.factorize(meta["_group"], sort=False)[0] if group_col else 0
    meta["_sample"] = meta.index
    meta = meta.sort_values(["_group_rank", "_label"], kind="stable")
    return meta


# --------------------------------------------------------------------------- #
# Panels
# --------------------------------------------------------------------------- #
def _style(opt: OverviewOptions) -> dict:
    return {
        "font.family": "DejaVu Sans", "font.size": opt.fs, "axes.labelsize": opt.fs,
        "axes.titlesize": opt.fs, "xtick.labelsize": opt.fs_small, "ytick.labelsize": opt.fs_small,
        "legend.fontsize": opt.fs_small, "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "svg.fonttype": "none", "pdf.fonttype": 42,
    }


def _clean(ax, keep=("left", "bottom")):
    for side, spine in ax.spines.items():
        spine.set_visible(side in keep)


class _Abundance:
    """(A) taxa x samples heatmap with metadata rows and bar charts."""

    def __init__(self, abund, tax, metadata, opt, colors):
        self.opt = opt
        self.meta = sample_metadata(metadata, list(abund.columns), opt.sample_label_col,
                                    opt.sample_group_col, opt.group_order)
        samples = self.meta["_sample"].tolist()
        heat, self.mags_per_taxon, mags_per_sample = abundance_by_taxon(
            abund, tax, opt.abundance_rank, opt.max_taxa)
        self.heat = heat[samples]
        self.mags_per_sample = mags_per_sample.reindex(samples)
        self.colors = colors
        meta_cols = opt.meta_cols or []
        self.meta_rows = []
        for col in meta_cols:
            if col not in self.meta.columns:
                raise ValueError(f"Metadata column {col!r} not found; available: {list(self.meta.columns)}")
            self.meta_rows.append((col, to_bool(self.meta[col])))

    def size_cm(self):
        o = self.opt
        cell = o.text_h_cm(o.fs_small) * 0.92
        n_s, n_t = self.heat.shape[1], self.heat.shape[0]
        label_w = max(len(t) for t in self.heat.index) * o.fs_small * PT_CM * 0.62
        width = label_w + n_s * cell + 0.2 + 1.4
        height = 1.3 + len(self.meta_rows) * cell * 0.85 + n_t * cell
        return width, height

    def draw(self, fig, spec):
        o = self.opt
        n_t, n_s = self.heat.shape
        n_m = len(self.meta_rows)
        cell = o.text_h_cm(o.fs_small) * 0.92
        g = spec.subgridspec(3, 2, height_ratios=[1.3, max(n_m, 0.001) * 0.85, n_t], width_ratios=[n_s * cell, 1.4],
                             hspace=0.05, wspace=0.03)
        ax_top = fig.add_subplot(g[0, 0])
        ax_meta = fig.add_subplot(g[1, 0], sharex=ax_top)
        ax_h = fig.add_subplot(g[2, 0], sharex=ax_top)
        ax_r = fig.add_subplot(g[2, 1], sharey=ax_h)
        self._cbar_slot = g[1, 1] if n_m else g[0, 1]
        self.axes = [ax_top, ax_meta, ax_h, ax_r]

        cmap = LinearSegmentedColormap.from_list("abund", ["#FFFFFF", "#FCBBA1", "#EF3B2C", "#67000D"])
        bounds = [0, 1, 2, 5, 10, 20, 40, 60, 80, 100]
        norm = BoundaryNorm(bounds, 256)
        vals = self.heat.values
        ax_h.imshow(vals, aspect="auto", cmap=cmap, norm=norm, extent=(0, n_s, n_t, 0), interpolation="nearest")
        for i in range(n_t):
            for j in range(n_s):
                v = vals[i, j]
                txt = f"{v:.0f}" if v >= 1 else ("<1" if v > 0 else "")
                ax_h.text(j + 0.5, i + 0.5, txt, ha="center", va="center", fontsize=o.fs_small * 0.82,
                          color="white" if v >= 40 else DARK)
        ax_h.set_yticks(np.arange(n_t) + 0.5)
        ax_h.set_yticklabels(self.heat.index)
        for lab in ax_h.get_yticklabels():
            lab.set_color(self.colors.get(lab.get_text(), DARK))
        ax_h.set_xticks(np.arange(n_s) + 0.5)
        ax_h.set_xticklabels(self.meta["_label"], rotation=90)
        ax_h.tick_params(length=0, pad=2)
        _clean(ax_h, keep=())
        groups = self.meta["_group"].tolist()
        for x in range(1, n_s):
            if groups[x] != groups[x - 1]:
                for ax in (ax_h, ax_meta):
                    ax.axvline(x, color="white", lw=1.6)

        # metadata rows: filled squares for booleans, coloured squares for categories
        ax_meta.set_xlim(0, n_s)
        ax_meta.set_ylim(max(n_m, 1), 0)
        legends = []
        for r, (col, as_bool) in enumerate(self.meta_rows):
            if as_bool is not None:
                for j, on in enumerate(as_bool.tolist()):
                    ax_meta.add_patch(Rectangle((j + 0.17, r + 0.12), 0.66, 0.76, lw=0.5, edgecolor="#9A9A9A",
                                                facecolor=DARK if on else "white"))
            else:
                cats = list(dict.fromkeys(self.meta[col].astype(str)))
                cmap_c = dict(zip(cats, CATEGORICAL[3:] + CATEGORICAL[:3]))
                for j, v in enumerate(self.meta[col].astype(str)):
                    ax_meta.add_patch(Rectangle((j + 0.1, r + 0.12), 0.8, 0.76, lw=0, facecolor=cmap_c[v]))
                legends.append((col, cmap_c))
        ax_meta.set_yticks(np.arange(n_m) + 0.5)
        ax_meta.set_yticklabels(o.meta_labels or [c for c, _ in self.meta_rows])
        ax_meta.tick_params(length=0, pad=2, labelbottom=False, bottom=False)
        _clean(ax_meta, keep=())
        if n_m == 0:
            ax_meta.set_visible(False)
        for col, cmap_c in legends:
            ax_r.legend(handles=[Patch(color=c, label=k) for k, c in cmap_c.items()], title=col,
                        loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=o.fs_small)

        # MAGs detected per sample (top) and MAGs per taxon (right)
        ax_top.bar(np.arange(n_s) + 0.5, self.mags_per_sample.values, width=0.72, color="#7A7A7A")
        ax_top.set_ylabel("MAGs\ndetected", labelpad=2)
        ax_top.set_xlim(0, n_s)
        ax_top.tick_params(axis="x", bottom=False, labelbottom=False)
        ax_top.yaxis.set_major_locator(plt.MaxNLocator(3, integer=True))
        _clean(ax_top)

        ax_r.barh(np.arange(n_t) + 0.5, self.mags_per_taxon.values, height=0.72,
                  color=[self.colors.get(t, OTHER_COLOR) for t in self.heat.index])
        for y, v in zip(np.arange(n_t) + 0.5, self.mags_per_taxon.values):
            ax_r.text(v, y, f" {v}", va="center", fontsize=o.fs_small * 0.9)
        ax_r.set_xlim(0, self.mags_per_taxon.max() * 1.35)
        ax_r.set_xlabel("MAGs", labelpad=2)
        ax_r.xaxis.set_major_locator(plt.MaxNLocator(2, integer=True))
        ax_r.tick_params(axis="y", left=False, labelleft=False)
        _clean(ax_r, keep=("bottom",))

        # colour bar below the right bar chart
        self._cbar = (cmap, norm)
        return ax_top

    def finish(self, fig):
        o = self.opt
        slot = self._cbar_slot.get_position(fig)
        h = o.text_h_cm() * 0.3 * CM / fig.get_figheight()
        cax = fig.add_axes([slot.x0 + slot.width * 0.08, slot.y0 + slot.height * 0.35, slot.width * 0.9, h])
        cmap, norm = self._cbar
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal",
                          ticks=[0, 10, 100])
        cb.outline.set_linewidth(0.4)
        cax.set_title("Abundance\n(%)", fontsize=o.fs_small, pad=2, linespacing=0.95)
        cax.tick_params(labelsize=o.fs_small * 0.9, length=1.5, pad=1)
        self.axes.append(cax)


class _Quality:
    """(B) completeness vs contamination with marginal histograms."""

    def __init__(self, quality: pd.DataFrame, opt: OverviewOptions):
        self.opt = opt
        q = quality[["Completeness", "Contamination"]].apply(pd.to_numeric, errors="coerce").dropna()
        self.q = q
        self.cat = quality_category(q.Completeness, q.Contamination, opt.hq_completeness,
                                    opt.hq_contamination, opt.mq_completeness, opt.mq_contamination)

    def size_cm(self):
        return 6.0, 6.0

    def draw(self, fig, spec):
        o = self.opt
        q, cat = self.q, self.cat
        g = spec.subgridspec(2, 2, width_ratios=[4.2, 1.1], height_ratios=[1.1, 4.2], wspace=0.06, hspace=0.06)
        ax = fig.add_subplot(g[1, 0])
        ax_hx = fig.add_subplot(g[0, 0], sharex=ax)
        ax_hy = fig.add_subplot(g[1, 1], sharey=ax)
        self.axes = [ax, ax_hx, ax_hy]
        ymax = o.max_contamination_axis
        xmin = min(o.mq_completeness, float(np.floor(q.Completeness.min() / 5) * 5)) if len(q) else 50
        xmin = max(0.0, xmin)
        colors = {"HQ": HQ_COLOR, "MQ": MQ_COLOR, "LQ": LQ_COLOR}
        ms = 2.6 * (o.fs / 8.0)
        for k in ["LQ", "MQ", "HQ"]:
            sel = cat == k
            inside = sel & (q.Contamination.values <= ymax)
            above = sel & (q.Contamination.values > ymax)
            if inside.any():
                ax.plot(q.Completeness.values[inside], q.Contamination.values[inside], "o", ms=ms,
                        color=colors[k], mew=0, zorder=3)
            if above.any():
                ax.plot(q.Completeness.values[above], np.full(above.sum(), ymax), "^", ms=ms * 1.25,
                        color=colors[k], mew=0, zorder=3, clip_on=False)
        ax.axvline(o.hq_completeness, ls="--", lw=0.7, color=HQ_COLOR)
        if o.mq_completeness > xmin:
            ax.axvline(o.mq_completeness, ls="--", lw=0.7, color=MQ_COLOR)
        ax.axhline(o.hq_contamination, ls="--", lw=0.7, color="#AAAAAA")
        ax.set_xlim(xmin, 100.5)
        ax.set_ylim(0, ymax)
        ax.set_xlabel("Completeness (%)")
        ax.set_ylabel("Contamination (%)", labelpad=2)
        ax.grid(True, ls=":", lw=0.4, alpha=0.7)
        _clean(ax)

        n = max(len(q), 1)
        bx = np.arange(xmin, 100.0 + 2.5, 2.5)
        by = np.arange(0, ymax + 0.5, 0.5)
        cy = q.Contamination.clip(upper=ymax - 1e-6)
        bottom_x = np.zeros(len(bx) - 1)
        left_y = np.zeros(len(by) - 1)
        for k in ["LQ", "MQ", "HQ"]:
            hx, _ = np.histogram(q.Completeness[cat == k], bins=bx)
            hy, _ = np.histogram(cy[cat == k], bins=by)
            ax_hx.bar(bx[:-1], hx / n * 100, width=np.diff(bx), align="edge", bottom=bottom_x, color=colors[k],
                      lw=0.3, edgecolor="white")
            ax_hy.barh(by[:-1], hy / n * 100, height=np.diff(by), align="edge", left=left_y, color=colors[k],
                       lw=0.3, edgecolor="white")
            bottom_x += hx / n * 100
            left_y += hy / n * 100
        ax_hx.set_ylabel("% MAGs", labelpad=2)
        ax_hx.yaxis.set_major_locator(plt.MaxNLocator(2, integer=True, prune="lower"))
        ax_hx.tick_params(axis="x", labelbottom=False)
        _clean(ax_hx)
        ax_hy.set_xlabel("% MAGs", labelpad=2)
        ax_hy.xaxis.set_major_locator(plt.MaxNLocator(2, integer=True, prune="lower"))
        ax_hy.tick_params(axis="y", labelleft=False)
        _clean(ax_hy)

        labels = {"HQ": "High quality", "MQ": "Medium quality", "LQ": "Other"}
        handles = [Line2D([], [], marker="o", ls="", color=colors[k], ms=ms * 1.3,
                          label=f"{labels[k]} (n = {(cat == k).sum()})") for k in ["HQ", "MQ", "LQ"]
                   if (cat == k).any()]
        above = q.Contamination.values > ymax
        if above.any():
            vals = ", ".join(f"{v:.1f}%" for v in sorted(q.Contamination.values[above]))
            handles.append(Line2D([], [], marker="^", ls="", color="#555555", ms=ms * 1.4,
                                  label=f"Contamination > {ymax:g}%\n({vals})"))
        ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(-0.02, -0.2), frameon=False, ncol=1,
                  handletextpad=0.1, borderpad=0.0, borderaxespad=0.0, labelspacing=0.25, fontsize=o.fs_small)
        return ax_hx


class _Taxonomy:
    """(C) taxonomic flow between two ranks, bands coloured by ``color_rank``."""

    def __init__(self, tax: pd.DataFrame, opt: OverviewOptions, colors):
        self.opt = opt
        a, b = (RANKS.index(r) for r in opt.flow_ranks)
        self.ranks = RANKS[a:b + 1]
        self.columns, self.values, self.children = taxonomy_flow(tax, self.ranks)
        self.colors = colors
        self.color_level = RANKS.index(opt.color_rank) - a  # index into node tuple, may be < 0
        self.tax = tax
        o = opt
        self.unit = o.text_h_cm(o.fs_small) * 0.36  # cm per genome
        self.gap = o.text_h_cm(o.fs_small) * 0.42   # cm between nodes
        self.col_h = [sum(self.values[n] for n in col) * self.unit + (len(col) - 1) * self.gap
                      for col in self.columns]

    def size_cm(self):
        return self.opt.width_cm, max(self.col_h) + self.opt.text_h_cm() * 1.6

    def _label(self, node, k):
        name = node[-1]
        name = "unclassified" if name == "Unclassified" else name
        show = getattr(self, "show_counts", None)
        return f"{name} ({self.values[node]})" if show is None or show[k] else name

    def _color(self, node):
        if self.color_level >= 0 and len(node) > self.color_level:
            return self.colors.get(node[self.color_level], OTHER_COLOR)
        # colour rank is above the first flow rank: look the lineage up
        name = node[0]
        col = self.ranks[0]
        hit = self.tax.loc[self.tax[col] == name, self.opt.color_rank]
        return self.colors.get(hit.iloc[0], OTHER_COLOR) if len(hit) else OTHER_COLOR

    def draw(self, fig, spec):
        o = self.opt
        ax = fig.add_subplot(spec)
        self.axes = [ax]
        width = o.width_cm
        H = max(self.col_h)
        char_w = o.fs_small * PT_CM * 0.6
        node_w = 0.13
        n_col = len(self.columns)

        def col_widths(counts_mid):
            self.show_counts = [True] + [counts_mid] * (n_col - 2) + [True] if n_col > 1 else [True]
            return [max(len(self._label(n, k)) for n in col) * char_w + 0.5 for k, col in enumerate(self.columns)]

        widths = col_widths(True)
        if sum(widths) + n_col * node_w > width:
            widths = col_widths(False)
        free = max(0.0, width - sum(widths) - n_col * node_w)
        xs, x = [], 0.0
        for k in range(n_col):
            xs.append(x)
            x += node_w + widths[k] + (free / (n_col - 1) if n_col > 1 else 0)
        pos = {}
        for k, col in enumerate(self.columns):
            y = (H - self.col_h[k]) / 2
            for node in col:
                h = self.values[node] * self.unit
                pos[node] = [y, y + h, y]
                y += h + self.gap

        def ribbon(x0, a0, b0, x1, a1, b1, color):
            xm = (x0 + x1) / 2
            verts = [(x0, a0), (xm, a0), (xm, a1), (x1, a1), (x1, b1), (xm, b1), (xm, b0), (x0, b0), (x0, a0)]
            codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4, MPath.LINETO, MPath.CURVE4,
                     MPath.CURVE4, MPath.CURVE4, MPath.CLOSEPOLY]
            ax.add_patch(PathPatch(MPath(verts, codes), facecolor=color, edgecolor="none", alpha=0.55, lw=0))

        for k in range(len(self.columns) - 1):
            for node in self.columns[k]:
                for child in self.children.get(node, []):
                    v = self.values[child] * self.unit
                    p = pos[node]
                    a0 = p[2]
                    p[2] += v
                    ct, cb, _ = pos[child]
                    ribbon(xs[k] + node_w, a0, a0 + v, xs[k + 1], ct, cb, self._color(child))
        th = o.text_h_cm(o.fs_small) * 0.95
        for k, col in enumerate(self.columns):
            centres = np.array([(pos[n][0] + pos[n][1]) / 2 for n in col])
            label_y = spread_labels(centres, th, -0.2, H + 0.2)
            for node, c, ly in zip(col, centres, label_y):
                t, b, _ = pos[node]
                ax.add_patch(Rectangle((xs[k], t), node_w, b - t, color=DARK, lw=0))
                last = k == len(self.columns) - 1
                x_text = xs[k] + node_w + (0.28 if abs(ly - c) > 1e-6 else 0.06)
                if abs(ly - c) > 1e-6:
                    ax.plot([xs[k] + node_w + 0.02, xs[k] + node_w + 0.24], [c, ly], lw=0.4, color="#777777")
                name = node[-1]
                ax.text(x_text, ly, self._label(node, k), va="center", ha="left", fontsize=o.fs_small,
                        fontstyle="italic" if self.ranks[k] in ("genus", "species") and name != "Unclassified"
                        else "normal",
                        bbox=None if last else dict(boxstyle="square,pad=0.04", fc="white", ec="none", alpha=0.65))
        for k, rank in enumerate(self.ranks):
            ax.text(xs[k] + node_w / 2, -o.text_h_cm() * 0.5, rank.capitalize(), ha="center", va="bottom",
                    fontsize=o.fs, fontweight="bold")
        ax.set_xlim(0, width)
        ax.set_ylim(H + 0.05, -o.text_h_cm() * 1.6)
        ax.axis("off")
        return ax


class _Function:
    """(D) largest clusters with taxonomy and KEGG module completeness."""

    def __init__(self, drep, gtdb, pathways, opt: OverviewOptions, colors):
        self.opt = opt
        self.colors = colors
        cluster_data, _, self.total_clusters = prepare_cluster_data(
            drep, gtdb, (opt.color_rank, "genus", "species"), opt.top_n * 2)
        pw_genomes = set(pathways["contig"].map(normalize_genome_id))
        self.n_without_kegg = int((~cluster_data.representative.isin(pw_genomes)).sum())
        cd = cluster_data[cluster_data.representative.isin(pw_genomes)]
        self.clusters = cd.sort_values("n_members", ascending=False, kind="stable").head(opt.top_n)
        self.clusters = self.clusters.reset_index(drop=True)
        reps = self.clusters.representative.tolist()
        hm, self.module_meta, self.blocks = prepare_pathway_matrix(pathways, reps, top_modules=opt.top_modules,
                                                                   mode=opt.module_mode)
        self.hm = hm.reindex(reps)
        self.n_mags = len(drep)

    def size_cm(self):
        o = self.opt
        row = o.text_h_cm(o.fs_small) * 0.95
        return o.width_cm, 1.0 + len(self.clusters) * row

    @staticmethod
    def species_label(row):
        s, g = row.get("species"), row.get("genus")
        if isinstance(s, str) and s and s != "Unclassified":
            return s
        return f"{g} sp." if isinstance(g, str) and g and g != "Unclassified" else "unclassified"

    def draw(self, fig, spec):
        o = self.opt
        labels = [self.species_label(r) for _, r in self.clusters.iterrows()]
        char_w = o.fs_small * PT_CM * 0.6
        txt_w = max(len(s) for s in labels) * char_w + 0.3
        bar_w = 2.2
        heat_w = o.width_cm - txt_w - bar_w
        self._heat_w = heat_w
        g = spec.subgridspec(2, 3, width_ratios=[bar_w, txt_w, heat_w], height_ratios=[1.0, 10], wspace=0.0,
                             hspace=0.02)
        ax_bar = fig.add_subplot(g[1, 0])
        ax_txt = fig.add_subplot(g[1, 1], sharey=ax_bar)
        ax_k = fig.add_subplot(g[1, 2], sharey=ax_bar)
        ax_cls = fig.add_subplot(g[0, 2], sharex=ax_k)
        self.axes = [ax_bar, ax_txt, ax_k, ax_cls]
        nr = len(self.clusters)
        yy = np.arange(nr) + 0.5
        n_members = self.clusters.n_members.values
        ax_bar.barh(yy, n_members, height=0.72,
                    color=[self.colors.get(p, OTHER_COLOR) for p in self.clusters[o.color_rank]])
        for y, v in zip(yy, n_members):
            ax_bar.text(v, y, f" {v}", va="center", fontsize=o.fs_small * 0.9)
        ax_bar.set_ylim(nr, 0)
        ax_bar.set_xlim(0, n_members.max() * 1.3)
        ax_bar.xaxis.set_major_locator(plt.MaxNLocator(3, integer=True))
        ax_bar.tick_params(axis="y", left=False, labelleft=False)
        ax_bar.set_xlabel("MAGs per cluster", labelpad=2)
        _clean(ax_bar, keep=("bottom",))

        ax_txt.set_xlim(0, 1)
        ax_txt.axis("off")
        for y, s in zip(yy, labels):
            ax_txt.text(0.04, y, s, va="center", fontsize=o.fs_small,
                        fontstyle="normal" if s == "unclassified" else "italic")
        ax_txt.text(0.04, -0.4, "Species", fontsize=o.fs_small, fontweight="bold", va="bottom")

        cmap = LinearSegmentedColormap.from_list("kegg", ["#FFFFFF", "#FEE0D2", "#FB6A4A", "#A50F15", "#67000D"])
        m = self.hm.shape[1]
        self._im = ax_k.imshow(self.hm.fillna(0).values, aspect="auto", cmap=cmap, vmin=0, vmax=100,
                               extent=(0, m, nr, 0), interpolation="nearest")
        for x in range(1, m):
            ax_k.axvline(x, color="white", lw=0.3)
        for y in range(1, nr):
            ax_k.axhline(y, color="white", lw=0.3)
        ax_k.set_xticks(np.arange(m) + 0.5)
        ax_k.set_xticklabels(self.module_meta.module_accession, rotation=90, fontsize=o.fs_small * 0.85)
        ax_k.tick_params(axis="y", left=False, labelleft=False)
        ax_k.tick_params(axis="x", length=0, pad=1.5)
        _clean(ax_k, keep=())

        ax_cls.set_ylim(0, 1)
        ax_cls.axis("off")
        # class labels: one line each; overlapping neighbours are moved up a level in finish()
        self._cls_texts = []
        for i, (cls, a, b) in enumerate(self.blocks):
            col = CLASS_COLORS[i % len(CLASS_COLORS)]
            ax_cls.add_patch(Rectangle((a + 0.08, 0.0), b - a + 0.84, 0.16, color=col, lw=0))
            if a > 0:
                ax_k.axvline(a, color="#555555", lw=0.7)
            label = CLASS_SHORT.get(cls, cls).replace("\n", " ")
            self._cls_texts.append(ax_cls.text((a + b + 1) / 2, 0.24, label, fontsize=o.fs_small * 0.9,
                                               ha="center", va="bottom", color=col))
        return ax_cls

    def finish(self, fig):
        o = self.opt
        r = fig.canvas.get_renderer()
        pad = o.fs_small * 0.5 * fig.dpi / 72
        ends: list[float] = []
        for t in self._cls_texts:
            bb = t.get_window_extent(r)
            level = next((i for i, e in enumerate(ends) if bb.x0 > e + pad), len(ends))
            if level == len(ends):
                ends.append(0.0)
            ends[level] = bb.x1
            t.set_y(0.24 + level * 0.42)
        ax_bar, ax_k = self.axes[0], self.axes[2]
        inv = fig.transFigure.inverted()
        y_lab = inv.transform((0, ax_k.get_tightbbox().y0))[1]
        bb = ax_k.get_position()
        h = o.text_h_cm() * 0.28 * CM / fig.get_figheight()
        cax = fig.add_axes([bb.x1 - bb.width * 0.3, y_lab - h - o.text_h_cm() * 1.3 * CM / fig.get_figheight(),
                            bb.width * 0.3, h])
        cb = fig.colorbar(self._im, cax=cax, orientation="horizontal", ticks=[0, 50, 100])
        cb.outline.set_linewidth(0.4)
        cax.set_title("Module completeness (%)", fontsize=o.fs_small, pad=2)
        cax.tick_params(labelsize=o.fs_small * 0.9, length=1.5, pad=1)
        shown = int(self.clusters.n_members.sum())
        y_bar = inv.transform((0, ax_bar.get_tightbbox().y0))[1]
        fig.text(ax_bar.get_position().x0, y_bar - 0.004,
                 f"{self.n_mags} MAGs in {self.total_clusters} clusters;\n"
                 f"shown: {shown} MAGs ({shown / max(self.n_mags, 1) * 100:.0f}%)",
                 fontsize=o.fs_small * 0.9, va="top", color="#555555")
        self.axes.append(cax)


# --------------------------------------------------------------------------- #
# Figure assembly
# --------------------------------------------------------------------------- #
def _content_bbox(fig, axes):
    """Union of tight bounding boxes of ``axes`` in cm (x0, y0, x1, y1), figure origin bottom-left."""
    r = fig.canvas.get_renderer()
    bbs = [a.get_tightbbox(r) for a in axes if a.get_visible()]
    bbs = [b for b in bbs if b is not None]
    x0 = min(b.x0 for b in bbs); y0 = min(b.y0 for b in bbs)
    x1 = max(b.x1 for b in bbs); y1 = max(b.y1 for b in bbs)
    k = 2.54 / fig.dpi
    return x0 * k, y0 * k, x1 * k, y1 * k


def overview_figure(output_prefix, opt: OverviewOptions, coverm=None, gtdb=None, quality=None, drep=None,
                    pathways=None, metadata=None, formats=("png",), dpi=300):
    """Draw the overview figure and save it as ``<output_prefix>.<fmt>``.

    Panels whose inputs are missing are skipped with a warning. Returns the list
    of written files.
    """
    if gtdb is None:
        raise ValueError("A GTDB-Tk summary (--gtdb) is required for every panel")
    tax = taxonomy_table(gtdb)
    want = [p for p in PANELS if p in opt.panels]
    need = {"abundance": coverm, "quality": quality, "taxonomy": True, "function": drep}
    if pathways is None and "function" in want:
        need["function"] = None
    panels = [p for p in want if need[p] is not None]
    for p in want:
        if p not in panels:
            print(f"[WARN] Skipping panel '{p}': required input missing")
    if not panels:
        raise ValueError("No panel can be drawn with the given inputs")

    # colours of the colour rank, ordered by abundance if CoverM is available
    if coverm is not None:
        abund = coverm_relative_abundance(coverm)
        totals = abund.groupby(tax.reindex(abund.index)[opt.color_rank].fillna("Unclassified")).sum().sum(axis=1)
        taxa_order = totals.sort_values(ascending=False).index.tolist()
    else:
        abund = None
        taxa_order = tax[opt.color_rank].value_counts().index.tolist()
    colors = taxon_colors(taxa_order + ["Other", "Unclassified"])

    with plt.rc_context(_style(opt)):
        objs = {}
        if "abundance" in panels:
            objs["abundance"] = _Abundance(abund, tax, metadata, opt, colors)
        if "quality" in panels:
            objs["quality"] = _Quality(quality, opt)
        if "taxonomy" in panels:
            objs["taxonomy"] = _Taxonomy(tax, opt, colors)
        if "function" in panels:
            objs["function"] = _Function(drep, gtdb, pathways, opt, colors)

        # rows: [abundance, quality] side by side, then taxonomy, then function
        rows = []
        top = [p for p in ("abundance", "quality") if p in objs]
        if top:
            rows.append(top)
        rows += [[p] for p in ("taxonomy", "function") if p in objs]
        heights = [max(objs[p].size_cm()[1] for p in row) for row in rows]

        gaps = [1.0] * (len(rows) - 1)  # cm between rows, refined after measuring
        margins = [0.3, 0.3, 0.6, 0.2]  # left, right, top, bottom in cm (refined)
        for _ in range(3):
            plt.close("all")
            fig, row_axes = _build(objs, rows, heights, gaps, margins, opt)
            fig.canvas.draw()
            # vertical: overhang of each row beyond its slot, then set gaps to a fixed clearance
            clear = opt.text_h_cm() * 0.9
            ext = [_content_bbox(fig, ax) for ax in row_axes]
            W, H = fig.get_figwidth() * 2.54, fig.get_figheight() * 2.54
            slots = _slots(heights, gaps, margins)
            new_gaps = []
            for i in range(len(rows) - 1):
                below_overhang = slots[i][0] - ext[i][1]       # row i content below its slot bottom
                above_overhang = ext[i + 1][3] - slots[i + 1][1]  # row i+1 content above its slot top
                new_gaps.append(max(0.1, below_overhang + above_overhang + clear))
            x0 = min(e[0] for e in ext); x1 = max(e[2] for e in ext)
            y0 = min(e[1] for e in ext); y1 = max(e[3] for e in ext)
            new_margins = [margins[0] - x0 + 0.05, margins[1] - (W - x1) + 0.05,
                           margins[2] - (H - y1) + 0.05 + opt.text_h_cm(opt.fs_letter),
                           margins[3] - y0 + 0.05]
            converged = np.allclose(new_gaps, gaps, atol=0.02) and np.allclose(new_margins, margins, atol=0.02)
            gaps, margins = new_gaps, new_margins
            if converged:
                break
        written = []
        for fmt in formats:
            path = f"{output_prefix}.{fmt}"
            fig.savefig(path, dpi=dpi)
            written.append(path)
            print(f"[INFO] Saved {path} ({fig.get_figwidth() * 2.54:.1f} x {fig.get_figheight() * 2.54:.1f} cm)")
        plt.close(fig)
    return written


def _slots(heights, gaps, margins):
    """Bottom/top (cm, from the figure bottom) of every row slot."""
    total = margins[2] + sum(heights) + sum(gaps) + margins[3]
    slots, y = [], total - margins[2]
    for i, h in enumerate(heights):
        slots.append((y - h, y))
        y -= h + (gaps[i] if i < len(gaps) else 0)
    return slots


def _build(objs, rows, heights, gaps, margins, opt):
    W = opt.width_cm
    total_h = margins[2] + sum(heights) + sum(gaps) + margins[3]
    fig = plt.figure(figsize=(W * CM, total_h * CM))
    ratios = []
    for i, h in enumerate(heights):
        ratios.append(h)
        if i < len(gaps):
            ratios.append(gaps[i])
    gs = fig.add_gridspec(len(ratios), 1, height_ratios=ratios, hspace=0.0,
                          left=margins[0] / W, right=1 - margins[1] / W,
                          top=1 - margins[2] / total_h, bottom=margins[3] / total_h)
    row_axes = []
    inner_w = W - margins[0] - margins[1]
    for i, row in enumerate(rows):
        spec = gs[2 * i]
        axes = []
        if len(row) == 2:
            wa, _ = objs[row[0]].size_cm()
            wb = inner_w - wa - 1.3
            sub = spec.subgridspec(1, 3, width_ratios=[wa, 1.3, max(wb, 4.0)], wspace=0.0)
            specs = [sub[0], sub[2]]
        else:
            specs = [spec]
        for p, s in zip(row, specs):
            obj = objs[p]
            obj.draw(fig, s)
            axes += obj.axes
        row_axes.append(axes)
    fig.canvas.draw()
    for i, row in enumerate(rows):
        for p in row:
            if hasattr(objs[p], "finish"):
                objs[p].finish(fig)
                row_axes[i] = list(row_axes[i]) + [objs[p].axes[-1]]
    # panel letters at the top-left of each panel's content
    if len(objs) > 1:
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        inv = fig.transFigure.inverted()
        for i, row in enumerate(rows):
            for p in row:
                axes = [a for a in objs[p].axes if a.get_visible()]
                bb = [a.get_tightbbox(r) for a in axes]
                x0 = min(b.x0 for b in bb)
                y1 = max(b.y1 for b in bb)
                x, y = inv.transform((x0, y1))
                if p == "abundance" or len(row) == 1:
                    x = margins[0] / W * 0 + 0.0
                fig.text(max(x - 0.004, 0.0), y + 0.002, PANEL_LETTERS[p], fontsize=opt.fs_letter,
                         fontweight="bold", va="bottom", ha="left")
    return fig, row_axes
