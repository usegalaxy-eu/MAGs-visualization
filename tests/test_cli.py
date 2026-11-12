from __future__ import annotations
import subprocess
import sys
from pathlib import Path
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
TEST_DATA = REPO_ROOT / "test-data"


def run_cli(args: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    cmd = [sys.executable, "-m", "mags_visualization.main"] + args
    return subprocess.run(
        cmd,
        cwd=str(cwd or REPO_ROOT),
        capture_output=True,
        text=True,
    )


@pytest.mark.parametrize(
    "args",
    [
        ["--help"],
        ["sample-heatmap", "--help"],
        ["comp-conta", "--help"],
        ["taxa-sankey", "--help"],
        ["drep-cluster-annot", "--help"],
        ["drep-cluster-func", "--help"],
        ["pathway-module-heatmap", "--help"],
        ["overview-figure", "--help"],
    ],
)
def test_help_commands_work(args: list[str]) -> None:
    result = run_cli(args)
    assert result.returncode == 0, result.stderr


def test_sample_heatmap_smoke(tmp_path: Path) -> None:
    outdir = tmp_path / "sample-heatmap"
    result = run_cli(
        [
            "sample-heatmap",
            "--coverm", str(TEST_DATA / "coverm.tsv"),
            "--gtdb", str(TEST_DATA / "gtdb.tsv"),
            "--metadata", str(TEST_DATA / "metadata.tsv"),
            "--meta_cols",
            "Infection by Nosema ceranae",
            "Chronic exposure to neonicotinoid",
            "Treatment with probiotic",
            "--tax_level", "phylum",
            "--output", str(outdir),
        ]
    )
    assert result.returncode == 0, result.stderr
    output_file = outdir / "heatmap_with_bars_phylum.png"
    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_comp_conta_smoke(tmp_path: Path) -> None:
    outdir = tmp_path / "comp-conta"
    result = run_cli(
        [
            "comp-conta",
            "--checkm", str(TEST_DATA / "checkm.tsv"),
            "--checkm2", str(TEST_DATA / "checkm2.tsv"),
            "--mode", "tax",
            "--gtdb", str(TEST_DATA / "gtdb.tsv"),
            "--tax_level", "phylum",
            "--output", str(outdir),
        ]
    )
    assert result.returncode == 0, result.stderr
    output_file = outdir / "comp_conta_by_rank_marginals_CheckM2.png"
    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_taxa_sankey_smoke(tmp_path: Path) -> None:
    outdir = tmp_path / "taxa-sankey"
    result = run_cli(
        [
            "taxa-sankey",
            "--gtdb", str(TEST_DATA / "gtdb.tsv"),
            "--output", str(outdir),
        ]
    )
    assert result.returncode == 0, result.stderr
    output_file = outdir / "sankey_plot.html"
    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_drep_cluster_smoke(tmp_path: Path) -> None:
    outdir = tmp_path / "drep-cluster-annot"
    result = run_cli(
        [
            "drep-cluster-annot",
            "--drep", str(TEST_DATA / "drep.csv"),
            "--gtdb", str(TEST_DATA / "gtdb.tsv"),
            "--checkm2", str(TEST_DATA / "checkm2.tsv"),
            "--top_n", "30",
            "--output", str(outdir),
        ]
    )
    assert result.returncode == 0, result.stderr
    output_file = outdir / "drep_cluster_top30_phylum-genus.png"
    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_functional_annotation_smoke(tmp_path: Path) -> None:
    outdir = tmp_path / "drep-cluster-func"
    result = run_cli(
        [
            "drep-cluster-func",
            "--drep", str(TEST_DATA / "drep.csv"),
            "--gtdb", str(TEST_DATA / "gtdb.tsv"),
            "--pathways", str(TEST_DATA / "kegg_pathway_completeness.tsv"),
            "--top_n", "30",
            "--output", str(outdir),
        ]
    )
    assert result.returncode == 0, result.stderr
    output_file = outdir / "drep_cluster_functional_core_top30.png"
    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_pathway_module_heatmap_smoke(tmp_path: Path) -> None:
    outdir = tmp_path / "pathway-module-heatmap"
    result = run_cli(
        [
            "pathway-module-heatmap",
            "--drep", str(TEST_DATA / "drep.csv"),
            "--gtdb", str(TEST_DATA / "gtdb.tsv"),
            "--pathways", str(TEST_DATA / "kegg_pathway_completeness.tsv"),
            "--format", "png",
            "-o", str(outdir),
        ]
    )
    assert result.returncode == 0, result.stderr + result.stdout
    assert any(outdir.glob("*.png"))


def test_taxa_sankey_abundance_and_ordered_smoke(tmp_path: Path) -> None:
    outdir = tmp_path / "taxa-sankey"
    result = run_cli(
        [
            "taxa-sankey",
            "--gtdb", str(TEST_DATA / "gtdb.tsv"),
            "--coverm", str(TEST_DATA / "coverm.tsv"),
            "--font_size", "14",
            "--order_links",
            "--format", "html",
            "-o", str(outdir),
        ]
    )
    assert result.returncode == 0, result.stderr + result.stdout
    # --order_links writes the crossing-minimised version instead of the plain one
    assert (outdir / "sankey_plot_ordered.html").exists()
    assert (outdir / "sankey_plot_rank_filtered.html").exists()


# ---- input validation ---- #

def test_sample_heatmap_requires_coverm() -> None:
    result = run_cli(
        [
            "sample-heatmap",
            "--gtdb", str(TEST_DATA / "gtdb.tsv"),
            "--output", "dummy",
        ]
    )
    assert result.returncode != 0


def test_sample_heatmap_requires_gtdb() -> None:
    result = run_cli(
        [
            "sample-heatmap",
            "--coverm", str(TEST_DATA / "coverm.tsv"),
            "--output", "dummy",
        ]
    )
    assert result.returncode != 0


def test_comp_conta_requires_gtdb_when_coloring_by_tax() -> None:
    result = run_cli(
        [
            "comp-conta",
            "--checkm", str(TEST_DATA / "checkm.tsv"),
            "--checkm2", str(TEST_DATA / "checkm2.tsv"),
            "--color_by", "tax",
            "--output", "dummy",
        ]
    )
    assert result.returncode != 0


def test_comp_conta_requires_metadata_when_coloring_by_meta() -> None:
    result = run_cli(
        [
            "comp-conta",
            "--checkm", str(TEST_DATA / "checkm.tsv"),
            "--checkm2", str(TEST_DATA / "checkm2.tsv"),
            "--color_by", "meta",
            "--meta_col", "Infection by Nosema ceranae",
            "--output", "dummy",
        ]
    )
    assert result.returncode != 0


def test_comp_conta_requires_meta_col_when_coloring_by_meta() -> None:
    result = run_cli(
        [
            "comp-conta",
            "--checkm", str(TEST_DATA / "checkm.tsv"),
            "--checkm2", str(TEST_DATA / "checkm2.tsv"),
            "--color_by", "meta",
            "--metadata", str(TEST_DATA / "metadata.tsv"),
            "--output", "dummy",
        ]
    )
    assert result.returncode != 0


def test_drep_cluster_requires_drep() -> None:
    result = run_cli(
        [
            "drep-cluster-annot",
            "--gtdb", str(TEST_DATA / "gtdb.tsv"),
            "--checkm2", str(TEST_DATA / "checkm2.tsv"),
            "--output", "dummy",
        ]
    )
    assert result.returncode != 0


def test_taxa_sankey_requires_gtdb() -> None:
    result = run_cli(
        [
            "taxa-sankey",
            "--output", "dummy",
        ]
    )
    assert result.returncode != 0


def test_functional_annotation_requires_kegg_pathway_completeness() -> None:
    result = run_cli(
        [
            "drep-cluster-func",
            "--drep", str(TEST_DATA / "drep.csv"),
            "--gtdb", str(TEST_DATA / "gtdb.tsv"),
            "--output", "dummy",
        ]
    )
    assert result.returncode != 0


def test_functional_annotation_requires_drep() -> None:
    result = run_cli(
        [
            "drep-cluster-func",
            "--gtdb", str(TEST_DATA / "gtdb.tsv"),
            "--kegg_pathway_completeness", str(TEST_DATA / "kegg_pathway_completeness.tsv"),
            "--output", "dummy",
        ]
    )
    assert result.returncode != 0
