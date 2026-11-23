
from __future__ import annotations
import numpy as np
import pandas as pd
import pytest

from mags_visualization.id_normalizer import format_genome_display, normalize_genome_id, normalize_genome_series
from mags_visualization.drep_cluster_func import extract_rank, format_genome_display, simplify_pathway_class
from mags_visualization.heatmap import clean_sample_from_coverm_col


def test_simplify_pathway_class_second_level() -> None:
    assert simplify_pathway_class("Metabolism; Carbohydrate metabolism") == "Carbohydrate metabolism"


def test_simplify_pathway_class_single_level() -> None:
    assert simplify_pathway_class("Metabolism") == "Metabolism"


def test_simplify_pathway_class_empty_string() -> None:
    assert simplify_pathway_class("") == "Unclassified"


def test_simplify_pathway_class_nan() -> None:
    assert simplify_pathway_class(np.nan) == "Unclassified"


def test_extract_rank_phylum() -> None:
    tax = "d__Bacteria; p__Pseudomonadota; c__Gammaproteobacteria"
    assert extract_rank(tax, "phylum") == "Pseudomonadota"


def test_extract_rank_missing_rank_unclassified() -> None:
    tax = "d__Bacteria; p__Pseudomonadota"
    assert extract_rank(tax, "species") == "Unclassified"


def test_extract_rank_nan_unclassified() -> None:
    assert extract_rank(np.nan, "genus") == "Unclassified"


def test_extract_rank_no_rank() -> None:
    with pytest.raises(ValueError):
        extract_rank("d__Bacteria; p__Firmicutes", "kingdom")


def test_extract_rank_empty_name_unclassified() -> None:
    tax = "d__Bacteria; g__"
    assert extract_rank(tax, "genus") == "Unclassified"


def test_normalize_genome_id_removes_extension() -> None:
    assert normalize_genome_id("SRR24759597_bin_73064.fasta") == "srr24759597_bin_73064_fasta"


def test_normalize_genome_id_removes_count_suffix() -> None:
    assert normalize_genome_id("ERR9966616_bin_10_fasta_Count") == "err9966616_bin_10_fasta"


def test_normalize_genome_id_paths() -> None:
    assert normalize_genome_id("/tmp/data/group_2_bin_200650.fna") == "group_2_bin_200650_fasta"


def test_normalize_genome_id_appends_fasta_if_missing() -> None:
    assert normalize_genome_id("group_2_bin_200650") == "group_2_bin_200650_fasta"


def test_normalize_genome_series() -> None:
    series = pd.Series(["SRR1_bin_1.fasta", "group_2_bin_200650"])
    result = normalize_genome_series(series)
    assert result.tolist() == ["srr1_bin_1_fasta", "group_2_bin_200650_fasta"]


def test_format_genome_display_srr() -> None:
    assert format_genome_display("srr24759597_bin_73064_fasta") == "SRR24759597_bin_73064"


def test_format_genome_display_group() -> None:
    assert format_genome_display("group_2_bin_200650_fasta") == "group_2_bin_200650"


def test_clean_sample_from_coverm_col_extracts_srr() -> None:
    col = "SRR12345678 something else"
    assert clean_sample_from_coverm_col(col) == "SRR12345678"


def test_clean_sample_from_coverm_col_removes_fastq() -> None:
    col = "sampleA.fastq extra"
    assert clean_sample_from_coverm_col(col) == "sampleA"

def test_gtdb_table_finds_user_genome_in_any_column(tmp_path):
    from mags_visualization.main import gtdb_table, read_table

    path = tmp_path / "gtdb.tsv"
    path.write_text("classification\tuser_genome\tclosest_genome_ani\n"
                    "d__Bacteria;p__Bacillota\tbin_1.fasta\t98.1\n")
    df = gtdb_table(read_table(str(path), index_col=None, prefer_tsv=True))
    assert df.index.name == "user_genome"
    assert list(df.index) == ["bin_1.fasta"]
    assert df.loc["bin_1.fasta", "classification"] == "d__Bacteria;p__Bacillota"
