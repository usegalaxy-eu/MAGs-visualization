"""Tests of the overview figure: data preparation, layout and the CLI subcommand."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mags_visualization.id_normalizer import normalize_genome_id
from mags_visualization.overview import (
    OverviewOptions,
    abundance_by_taxon,
    coverm_relative_abundance,
    overview_figure,
    quality_category,
    sample_metadata,
    spread_labels,
    taxon_colors,
    taxonomy_flow,
    taxonomy_table,
    to_bool,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
TEST_DATA = REPO_ROOT / "test-data"
META_COLS = ["Infection by Nosema ceranae", "Chronic exposure to neonicotinoid", "Treatment with probiotic"]


@pytest.fixture(scope="module")
def data():
    return {
        "gtdb": pd.read_csv(TEST_DATA / "gtdb.tsv", sep="\t"),
        "coverm": pd.read_csv(TEST_DATA / "coverm.tsv", sep="\t"),
        "quality": pd.read_csv(TEST_DATA / "checkm2.tsv", sep="\t"),
        "drep": pd.read_csv(TEST_DATA / "drep.csv"),
        "pathways": pd.read_csv(TEST_DATA / "kegg_pathway_completeness.tsv", sep="\t"),
        "metadata": pd.read_csv(TEST_DATA / "metadata.tsv", sep="\t"),
    }


# --------------------------------------------------------------------------- data preparation
def test_taxonomy_table_fills_missing_ranks():
    gtdb = pd.DataFrame({
        "user_genome": ["a.fasta", "b"],
        "classification": ["d__Bacteria;p__Bacillota;c__Bacilli;o__;f__;g__;s__",
                           "d__Bacteria;p__Pseudomonadota;c__Gamma;o__Ent;f__Entb;g__Gill;s__Gill apis"],
    })
    tax = taxonomy_table(gtdb)
    assert list(tax.columns) == ["domain", "phylum", "class", "order", "family", "genus", "species"]
    assert tax.loc[normalize_genome_id("a.fasta"), "order"] == "Unclassified"
    assert tax.loc[normalize_genome_id("b"), "species"] == "Gill apis"


def test_coverm_drops_unmapped_and_suffixes():
    cov = pd.DataFrame({"Genome": ["bin_1", "unmapped"], "S1 Relative Abundance (%)": [40.0, 60.0]})
    out = coverm_relative_abundance(cov)
    assert list(out.index) == [normalize_genome_id("bin_1")]
    assert list(out.columns) == ["S1"]


def test_abundance_by_taxon_pools_other():
    abund = pd.DataFrame({"S1": [10, 5, 1, 0.5], "S2": [0, 7, 2, 1]}, index=list("abcd"))
    tax = pd.DataFrame({"phylum": ["P1", "P2", "P3", "Unclassified"]}, index=list("abcd"))
    heat, per_taxon, per_sample = abundance_by_taxon(abund, tax, "phylum", max_taxa=2)
    assert list(heat.index) == ["P2", "P1", "Other", "Unclassified"]
    assert heat.loc["Other", "S2"] == 2
    assert per_taxon["Other"] == 1
    assert per_sample.to_dict() == {"S1": 4, "S2": 3}


def test_quality_category_thresholds():
    cat = quality_category([95, 95, 75, 60], [1, 6, 9, 1])
    assert list(cat) == ["HQ", "MQ", "MQ", "LQ"]
    cat = quality_category([75], [12], mq_contamination=10)
    assert list(cat) == ["LQ"]


def test_taxonomy_flow_is_a_tree_in_drawing_order():
    tax = pd.DataFrame({
        "phylum": ["A", "A", "A", "B"],
        "class": ["a1", "a1", "a2", "b1"],
        "order": ["x", "y", "x", "z"],
    })
    columns, values, children = taxonomy_flow(tax, ["phylum", "class", "order"])
    assert columns[0] == [("A",), ("B",)]
    assert values[("A",)] == 3
    # 'x' under a1 and 'x' under a2 are different nodes
    assert ("A", "a1", "x") in values and ("A", "a2", "x") in values
    # children follow their parent: all of A's grandchildren come before B's
    order_col = columns[2]
    assert order_col.index(("B", "b1", "z")) == len(order_col) - 1
    assert children[("A",)] == [("A", "a1"), ("A", "a2")]


def test_spread_labels_keeps_order_gap_and_bounds():
    y = np.array([1.0, 1.05, 1.1, 5.0])
    out = spread_labels(y, 0.3, 0.0, 6.0)
    assert np.all(np.diff(out) >= 0.3 - 1e-9)
    assert out[3] == pytest.approx(5.0)
    assert out[:3].mean() == pytest.approx(y[:3].mean())
    edge = spread_labels(np.array([0.0, 0.0]), 0.5, 0.0, 10.0)
    assert edge.min() >= 0.0


def test_taxon_colors_grey_for_other():
    c = taxon_colors(["P1", "Other", "P2", "Unclassified"])
    assert c["Other"] == c["Unclassified"] != c["P1"] != c["P2"]


def test_to_bool_and_sample_metadata(data):
    assert to_bool(pd.Series(["TRUE", "false", "True"])).tolist() == [True, False, True]
    assert to_bool(pd.Series(["a", "b"])) is None
    samples = coverm_relative_abundance(data["coverm"]).columns.tolist()
    meta = sample_metadata(data["metadata"], samples, label_col="Label", group_col="Condition",
                           group_order=["C", "P", "T", "TP", "N", "NP", "NT"])
    assert meta["_label"].tolist()[:3] == ["C1", "C2", "C3"]
    assert meta["_label"].tolist()[-1] == "NT3"
    with pytest.raises(ValueError, match="not found"):
        sample_metadata(data["metadata"], samples, label_col="nope")


def test_options_validation():
    with pytest.raises(ValueError, match="Unknown panel"):
        OverviewOptions(panels=["nope"])
    with pytest.raises(ValueError, match="higher to a lower rank"):
        OverviewOptions(flow_ranks=("genus", "phylum"))
    with pytest.raises(ValueError, match="meta_labels"):
        OverviewOptions(meta_cols=["a", "b"], meta_labels=["x"])


# --------------------------------------------------------------------------- figure
def _figure_size_cm(png: Path):
    from PIL import Image

    with Image.open(png) as im:
        w, h = im.size
        dpi = im.info.get("dpi", (300, 300))[0]
    return w / dpi * 2.54, h / dpi * 2.54


def test_full_overview_figure(tmp_path, data):
    opt = OverviewOptions(meta_cols=META_COLS, meta_labels=["N. ceranae", "Neonicotinoid", "Probiotic"],
                          sample_label_col="Label", sample_group_col="Condition")
    files = overview_figure(tmp_path / "fig", opt, coverm=data["coverm"], gtdb=data["gtdb"],
                            quality=data["quality"], drep=data["drep"], pathways=data["pathways"],
                            metadata=data["metadata"], formats=("png", "svg"))
    assert [Path(f).suffix for f in files] == [".png", ".svg"]
    width, height = _figure_size_cm(tmp_path / "fig.png")
    assert width == pytest.approx(17.0, abs=0.1)
    assert 15 < height < 40
    svg = (tmp_path / "fig.svg").read_text()
    # editable text in SVG, taxon names present
    assert "Pseudomonadota" in svg and "<text" in svg


@pytest.mark.parametrize("panel", ["abundance", "quality", "taxonomy", "function"])
def test_single_panels(tmp_path, data, panel):
    opt = OverviewOptions(panels=[panel], width_cm=12, font_size=7)
    overview_figure(tmp_path / panel, opt, coverm=data["coverm"], gtdb=data["gtdb"], quality=data["quality"],
                    drep=data["drep"], pathways=data["pathways"])
    width, _ = _figure_size_cm(tmp_path / f"{panel}.png")
    assert width == pytest.approx(12.0, abs=0.1)


def test_missing_inputs_skip_panels(tmp_path, data, capsys):
    overview_figure(tmp_path / "fig", OverviewOptions(), gtdb=data["gtdb"], quality=data["quality"])
    out = capsys.readouterr().out
    assert "Skipping panel 'abundance'" in out and "Skipping panel 'function'" in out
    assert (tmp_path / "fig.png").exists()


def test_no_drawable_panel_raises(data):
    with pytest.raises(ValueError, match="No panel"):
        overview_figure("x", OverviewOptions(panels=["function"]), gtdb=data["gtdb"])


def test_larger_font_makes_taller_figure(tmp_path, data):
    sizes = []
    for fs in (7, 10):
        overview_figure(tmp_path / f"f{fs}", OverviewOptions(panels=["taxonomy"], font_size=fs), gtdb=data["gtdb"])
        sizes.append(_figure_size_cm(tmp_path / f"f{fs}.png")[1])
    assert sizes[1] > sizes[0]


# --------------------------------------------------------------------------- CLI
def run_cli(args):
    return subprocess.run([sys.executable, "-m", "mags_visualization.main", *args], cwd=REPO_ROOT,
                          capture_output=True, text=True)


def test_cli_overview_figure(tmp_path):
    out = tmp_path / "out"
    res = run_cli([
        "overview-figure", "--gtdb", str(TEST_DATA / "gtdb.tsv"), "--coverm", str(TEST_DATA / "coverm.tsv"),
        "--checkm2", str(TEST_DATA / "checkm2.tsv"), "--drep", str(TEST_DATA / "drep.csv"),
        "--pathways", str(TEST_DATA / "kegg_pathway_completeness.tsv"),
        "--metadata", str(TEST_DATA / "metadata.tsv"), "--meta_cols", "7", "8", "9",
        "--sample_label_col", "10", "--sample_group_col", "5", "--top_n", "15",
        "--format", "png", "pdf", "--prefix", "bee", "-o", str(out),
    ])
    assert res.returncode == 0, res.stderr + res.stdout
    assert (out / "bee.png").exists() and (out / "bee.pdf").exists()


def test_cli_overview_uses_checkm_without_checkm2(tmp_path):
    res = run_cli(["overview-figure", "--gtdb", str(TEST_DATA / "gtdb.tsv"), "--checkm", str(TEST_DATA / "checkm.tsv"),
                   "--panels", "quality", "-o", str(tmp_path)])
    assert res.returncode == 0, res.stderr + res.stdout
    assert (tmp_path / "overview_figure.png").exists()


def test_cli_overview_bad_metadata_column(tmp_path):
    res = run_cli(["overview-figure", "--gtdb", str(TEST_DATA / "gtdb.tsv"), "--coverm", str(TEST_DATA / "coverm.tsv"),
                   "--metadata", str(TEST_DATA / "metadata.tsv"), "--sample_label_col", "nope",
                   "--panels", "abundance", "-o", str(tmp_path)])
    assert res.returncode != 0
    assert "not found" in res.stdout + res.stderr


# --------------------------------------------------------------------------- --style overview of single plots
@pytest.mark.parametrize("args, expected", [
    (["sample-heatmap", "--coverm", "coverm.tsv", "--gtdb", "gtdb.tsv", "--metadata", "metadata.tsv",
      "--meta_cols", "7", "8", "9", "--sample_label_col", "Label", "--sample_group_col", "Condition"],
     ["heatmap_overview_phylum.png"]),
    (["comp-conta", "--tools", "both", "--checkm", "checkm.tsv", "--checkm2", "checkm2.tsv", "--mq_contamination", "10"],
     ["comp_conta_overview_checkm.png", "comp_conta_overview_checkm2.png"]),
    (["taxa-sankey", "--gtdb", "gtdb.tsv", "--flow_ranks", "class", "species", "--format", "svg"],
     ["sankey_plot_overview.svg"]),
    (["drep-cluster-func", "--drep", "drep.csv", "--gtdb", "gtdb.tsv", "--pathways", "kegg_pathway_completeness.tsv",
      "--mode", "both", "--top_n", "10"],
     ["drep_cluster_functional_core_overview_top10.png", "drep_cluster_functional_difference_overview_top10.png"]),
])
def test_single_plots_overview_style(tmp_path, args, expected):
    files = {"coverm.tsv", "gtdb.tsv", "metadata.tsv", "checkm.tsv", "checkm2.tsv", "drep.csv",
             "kegg_pathway_completeness.tsv"}
    args = [str(TEST_DATA / a) if a in files else a for a in args]
    res = run_cli([*args, "--style", "overview", "--width_cm", "10", "-o", str(tmp_path)])
    assert res.returncode == 0, res.stderr + res.stdout
    for name in expected:
        assert (tmp_path / name).exists(), sorted(p.name for p in tmp_path.iterdir())
    png = [tmp_path / n for n in expected if n.endswith(".png")]
    for f in png:
        assert _figure_size_cm(f)[0] == pytest.approx(10.0, abs=0.1)


@pytest.mark.parametrize("args, message", [
    (["comp-conta", "--tools", "checkm2", "--checkm2", "checkm2.tsv", "--mode", "tax", "--gtdb", "gtdb.tsv"],
     "--mode quality only"),
    (["taxa-sankey", "--gtdb", "gtdb.tsv", "--format", "html"], "png, pdf or svg"),
])
def test_overview_style_rejects_unsupported_modes(tmp_path, args, message):
    args = [str(TEST_DATA / a) if a.endswith((".tsv", ".csv")) else a for a in args]
    res = run_cli([*args, "--style", "overview", "-o", str(tmp_path)])
    assert res.returncode != 0
    assert message in res.stdout + res.stderr


def test_classic_style_is_default(tmp_path):
    res = run_cli(["sample-heatmap", "--coverm", str(TEST_DATA / "coverm.tsv"), "--gtdb", str(TEST_DATA / "gtdb.tsv"),
                   "-o", str(tmp_path)])
    assert res.returncode == 0, res.stderr
    assert (tmp_path / "heatmap_with_bars_phylum.png").exists()


def test_quality_panel_without_gtdb(tmp_path, data):
    overview_figure(tmp_path / "q", OverviewOptions(panels=["quality"]), quality=data["quality"])
    assert (tmp_path / "q.png").exists()
    with pytest.raises(ValueError, match="GTDB"):
        overview_figure(tmp_path / "t", OverviewOptions(panels=["taxonomy"]))
