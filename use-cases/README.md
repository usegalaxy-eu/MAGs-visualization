# Use cases

Real MAG catalogues to try the plots on, each produced with the Galaxy MAG-generation and annotation
workflows. Every folder holds the inputs that `mags-visualization` reads (`gtdb.tsv`, `checkm2.tsv`,
`coverm.tsv`, `drep.csv`, `kegg_pathway_completeness.tsv`, `quast.tsv`, `bakta.tsv`, and `metadata.tsv`
where samples have metadata).

| Use case | Data | Samples | MAGs | Metadata columns |
|---|---|---|---|---|
| Bee gut microbiome | [`../test-data`](../test-data) (also the test data) | 21 | 161 | `Infection by Nosema ceranae`, `Chronic exposure to neonicotinoid`, `Treatment with probiotic` |
| Clouds | [`cloud-use-case/data`](cloud-use-case/data) | 14 | 230 | `Geographic origin of the air mass`, `Condition`, `Season` |
| Marine | [`marine-use-case/data`](marine-use-case/data) | 3 | 153 | – |
| Termite gut | [`termite-use-case/data`](termite-use-case/data) | 28 | 110 | `Species`, `Casts`, `Colony` |

Where available, `data/README.md` links the Galaxy history of the workflow run.

## Making the plots

Run from the repository root. `d` is the data folder of the use case, the commands are the same for all
of them (leave out `--metadata`/`--meta_cols` for the marine use case):

```bash
d=use-cases/termite-use-case/data
meta=(--metadata $d/metadata.tsv --meta_cols Species Casts Colony)

# page-width overview figure
mags-visualization overview-figure \
  --gtdb $d/gtdb.tsv --coverm $d/coverm.tsv --checkm2 $d/checkm2.tsv \
  --drep $d/drep.csv --pathways $d/kegg_pathway_completeness.tsv "${meta[@]}" \
  --format png pdf -o plots/termite

# the single plots
mags-visualization sample-heatmap --coverm $d/coverm.tsv --gtdb $d/gtdb.tsv "${meta[@]}" -o plots/termite
mags-visualization comp-conta --tools checkm2 --checkm2 $d/checkm2.tsv --mode quality -o plots/termite
mags-visualization taxa-sankey --gtdb $d/gtdb.tsv --format png -o plots/termite
mags-visualization drep-cluster-annot --drep $d/drep.csv --gtdb $d/gtdb.tsv \
  --checkm2 $d/checkm2.tsv --quast $d/quast.tsv --bakta $d/bakta.tsv -o plots/termite
mags-visualization drep-cluster-func --drep $d/drep.csv --gtdb $d/gtdb.tsv \
  --pathways $d/kegg_pathway_completeness.tsv -o plots/termite
mags-visualization pathway-module-heatmap --drep $d/drep.csv --gtdb $d/gtdb.tsv \
  --pathways $d/kegg_pathway_completeness.tsv -o plots/termite
```

Large catalogues give a tall taxonomy panel in the overview figure; stop it at a higher rank
(`--flow_ranks phylum family`) or show fewer clusters (`--top_n 20`). The command for the bee gut
figure of the FAIRyMAGs manuscript is the first example of the [gallery](../docs/plots.md).
