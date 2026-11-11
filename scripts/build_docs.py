#!/usr/bin/env python
"""Build the documentation: the CLI reference and the plot gallery.

    python scripts/build_docs.py            # docs/cli.md + docs/plots.md + docs/gallery/
    python scripts/build_docs.py --cli      # only docs/cli.md (fast, no plotting)
    python scripts/build_docs.py --check    # exit 1 if docs/cli.md is out of date

Every example is run through the installed CLI on the bundled ``test-data``
(the bee gut microbiome use case), so the gallery doubles as an end-to-end test
of every subcommand. Example commands are written with paths relative to the
repository root and can be copied as they are.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
GALLERY = DOCS / "gallery"
MAX_WIDTH_PX = 1600  # gallery images are downscaled to keep the repository small

META = ["Infection by Nosema ceranae", "Chronic exposure to neonicotinoid", "Treatment with probiotic"]
T = "test-data/"

# name, title, description, arguments (without -o)
EXAMPLES = [
    ("overview-figure", "Overview figure (all panels)",
     "Page-width figure combining relative abundance per sample, MAG quality, taxonomic composition and "
     "the functional profile of the largest species-level clusters. Drawn at its print size "
     "(`--width_cm`, `--font_size`), so the text is legible as printed.",
     ["overview-figure", "--gtdb", T + "gtdb.tsv", "--coverm", T + "coverm.tsv", "--checkm2", T + "checkm2.tsv",
      "--drep", T + "drep.csv", "--pathways", T + "kegg_pathway_completeness.tsv",
      "--metadata", T + "metadata.tsv", "--meta_cols", *META,
      "--meta_labels", "N. ceranae", "Neonicotinoid", "Probiotic",
      "--sample_label_col", "Label", "--sample_group_col", "Condition",
      "--group_order", "C", "P", "T", "TP", "N", "NP", "NT"]),
    ("overview-figure-single-panel", "Overview figure, single panel",
     "Any subset of panels can be drawn with `--panels`; here only the quality panel at half page width.",
     ["overview-figure", "--gtdb", T + "gtdb.tsv", "--checkm2", T + "checkm2.tsv", "--panels", "quality",
      "--width_cm", "8.5", "--font_size", "7"]),
    ("sample-heatmap", "Sample heatmap",
     "Relative abundance of each phylum per sample, with metadata bars, MAGs per taxon (top) and MAGs per "
     "sample (right).",
     ["sample-heatmap", "--coverm", T + "coverm.tsv", "--gtdb", T + "gtdb.tsv", "--metadata", T + "metadata.tsv",
      "--meta_cols", *META, "--tax_level", "phylum"]),
    ("comp-conta-quality", "Completeness vs contamination, coloured by quality",
     "CheckM and/or CheckM2 completeness against contamination with marginal histograms.",
     ["comp-conta", "--tools", "checkm2", "--checkm2", T + "checkm2.tsv", "--mode", "quality"]),
    ("comp-conta-tax", "Completeness vs contamination, coloured by taxonomy",
     "Points coloured by the GTDB rank given with `--tax_level`.",
     ["comp-conta", "--tools", "checkm2", "--checkm2", T + "checkm2.tsv", "--mode", "tax",
      "--gtdb", T + "gtdb.tsv", "--tax_level", "phylum"]),
    ("comp-conta-meta", "Completeness vs contamination, coloured by metadata",
     "Points coloured by a sample metadata column (`--meta_col`).",
     ["comp-conta", "--tools", "checkm2", "--checkm2", T + "checkm2.tsv", "--mode", "meta",
      "--metadata", T + "metadata.tsv", "--meta_col", "Condition"]),
    ("drep-cluster-annot", "dRep clusters with taxonomy and assembly statistics",
     "The largest species-level clusters with their size, taxonomy and CheckM2/QUAST/Bakta annotation.",
     ["drep-cluster-annot", "--drep", T + "drep.csv", "--gtdb", T + "gtdb.tsv", "--checkm2", T + "checkm2.tsv",
      "--quast", T + "quast.tsv", "--bakta", T + "bakta.tsv", "--top_n", "30"]),
    ("drep-cluster-func", "dRep clusters with KEGG module completeness",
     "The largest clusters with the completeness of the KEGG modules best covered by their representatives.",
     ["drep-cluster-func", "--drep", T + "drep.csv", "--gtdb", T + "gtdb.tsv",
      "--pathways", T + "kegg_pathway_completeness.tsv", "--top_n", "30"]),
    ("pathway-module-heatmap", "KEGG module heatmap",
     "KEGG module completeness of all MAGs, grouped by pathway class.",
     ["pathway-module-heatmap", "--drep", T + "drep.csv", "--gtdb", T + "gtdb.tsv",
      "--pathways", T + "kegg_pathway_completeness.tsv"]),
    ("taxa-sankey", "Taxonomy Sankey diagram",
     "Interactive (HTML) or static Sankey diagram of the GTDB taxonomy; `--order_links` minimises crossings "
     "and `--coverm` sizes the bands by abundance.",
     ["taxa-sankey", "--gtdb", T + "gtdb.tsv", "--order_links", "--font_size", "14", "--format", "png"]),
]


def cli_reference() -> str:
    sys.path.insert(0, str(ROOT))
    from mags_visualization.main import build_parser  # noqa: E402

    parser = build_parser()
    parser.prog = "mags-visualization"
    out = ["# Command-line reference", "",
           "<!-- generated by scripts/build_docs.py, do not edit by hand -->", "",
           "All subcommands are run as `mags-visualization <subcommand> [options]`. Defaults are shown in "
           "brackets. See [plots.md](plots.md) for an example of every plot.", "",
           "```text", parser.format_help().rstrip(), "```", ""]
    sub = next(a for a in parser._actions if isinstance(a, argparse._SubParsersAction))
    for name, p in sub.choices.items():
        p.prog = f"mags-visualization {name}"
        out += [f"## {name}", "", "```text", p.format_help().rstrip(), "```", ""]
    return "\n".join(out)


def run_example(name, args) -> list[Path]:
    out = GALLERY / name
    if out.exists():
        shutil.rmtree(out)
    cmd = [sys.executable, "-m", "mags_visualization.main", *args, "-o", str(out.relative_to(ROOT))]
    res = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    if res.returncode != 0:
        raise SystemExit(f"[ERROR] {name} failed:\n{res.stdout}\n{res.stderr}")
    # keep only the images, downscaled
    from PIL import Image

    images = []
    for f in sorted(out.rglob("*")):
        if f.suffix == ".png":
            with Image.open(f) as im:
                if im.width > MAX_WIDTH_PX:
                    im = im.resize((MAX_WIDTH_PX, round(im.height * MAX_WIDTH_PX / im.width)), Image.LANCZOS)
                im.save(f, optimize=True)
            images.append(f)
        elif f.is_file():
            f.unlink()
    print(f"[OK] {name}: {len(images)} image(s)")
    return images


def quote(a: str) -> str:
    return f'"{a}"' if (" " in a or not a) else a


def plots_page(results) -> str:
    out = ["# Plot gallery", "",
           "<!-- generated by scripts/build_docs.py, do not edit by hand -->", "",
           "Every plot of `mags-visualization`, made from the bundled `test-data` (bee gut microbiome, "
           "38 species-level representative MAGs). Run the commands from the repository root; options are "
           "described in [cli.md](cli.md).", ""]
    out += ["| Example | Subcommand |", "|---|---|"]
    out += [f"| [{title}](#{name}) | `{args[0]}` |" for name, title, _, args, _ in results]
    out.append("")
    for name, title, desc, args, images in results:
        cmd = "mags-visualization " + " \\\n  ".join(
            " ".join(quote(a) for a in grp) for grp in _group(args)) + f" \\\n  -o {name}"
        out += [f'<a id="{name}"></a>', f"## {title}", "", desc, "", "```bash", cmd, "```", ""]
        for img in images:
            out += [f"![{img.stem}]({img.relative_to(DOCS).as_posix()})", ""]
    return "\n".join(out)


def _group(args):
    """Split an argument list into one line per option for readable commands."""
    groups, cur = [], [args[0]]
    for a in args[1:]:
        if a.startswith("-"):
            groups.append(cur)
            cur = [a]
        else:
            cur.append(a)
    groups.append(cur)
    return groups


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cli", action="store_true", help="only write docs/cli.md")
    ap.add_argument("--check", action="store_true", help="fail if docs/cli.md is out of date")
    a = ap.parse_args()
    text = cli_reference()
    if a.check:
        current = (DOCS / "cli.md").read_text() if (DOCS / "cli.md").exists() else ""
        if current != text:
            print("[ERROR] docs/cli.md is out of date; run: python scripts/build_docs.py --cli")
            return 1
        print("[OK] docs/cli.md is up to date")
        return 0
    DOCS.mkdir(exist_ok=True)
    (DOCS / "cli.md").write_text(text)
    print("[OK] docs/cli.md")
    if a.cli:
        return 0
    results = [(n, t, d, args, run_example(n, args)) for n, t, d, args in EXAMPLES]
    (DOCS / "plots.md").write_text(plots_page(results))
    print("[OK] docs/plots.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
