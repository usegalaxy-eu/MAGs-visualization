# Example use case
All commands below assume that the command-line tool `mags-visualization` is installed and available in your environment.

Each plot is generated using a dedicated **subcommand**.

The working directory should be **MAGs-visualization/use-cases** 
so that the example paths resolve correctly.

If the command is not available, install the package first: ['README.md'](../README.md)

# Bash
## Bee Gut Microbiome

```bash
use_case_folder="bee-use-case/data"
```

### Overview figure (manuscript Figure 4)

```bash
mags-visualization overview-figure \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --coverm "$use_case_folder/coverm.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --drep "$use_case_folder/drep.csv" \
  --pathways "$use_case_folder/kegg_pathway_completeness.tsv" \
  --metadata "$use_case_folder/metadata.tsv" \
  --meta_cols "Infection by Nosema ceranae" "Chronic exposure to neonicotinoid" "Treatment with probiotic" \
  --meta_labels "N. ceranae" Neonicotinoid Probiotic \
  --sample_label_col Label --sample_group_col Condition --group_order C P T TP N NP NT \
  --format png svg pdf \
  -o bee-use-case/plots
```

### Sample heatmap

* Genus

```bash
mags-visualization sample-heatmap \
  --coverm "$use_case_folder/coverm.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --metadata "$use_case_folder/metadata.tsv" \
  --meta_cols "Infection by Nosema ceranae" "Chronic exposure to neonicotinoid" "Treatment with probiotic" \
  --tax_level genus \
  --top_bar_spacer -1.0 \
  --spacer_meta 2.5 \
  --output "$use_case_folder/plots" \
  --no_log \
  --top_bar_height 2.0
```

* Species

```bash
mags-visualization sample-heatmap \
  --coverm "$use_case_folder/coverm.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --metadata "$use_case_folder/metadata.tsv" \
  --meta_cols "Infection by Nosema ceranae" "Chronic exposure to neonicotinoid" "Treatment with probiotic" \
  --tax_level species \
  --top_bar_spacer -1.0 \
  --spacer_meta 2.5 \
  --output "$use_case_folder/plots" \
  --no_log \
  --top_bar_height 2.0
```

### dRep cluster plot with annotation (QUAST, BAKTA)

```bash
mags-visualization drep-cluster-annot \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --drep "$use_case_folder/drep.csv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --bakta "$use_case_folder/bakta.tsv" \
  --quast "$use_case_folder/quast.tsv" \
  --output out
```

### dRep cluster plot with functional annotation (KEGG pathway completeness)

```bash
mags-visualization drep-cluster-func \
  --drep "$use_case_folder/drep.csv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --pathways "$use_case_folder/kegg_pathway_completeness.tsv" \
  --top_modules 30 \
  --mode mean \
  --output out
  ```

To include all modules, use `--top_modules None`.

```bash
mags-visualization drep-cluster-func \
  --drep "$use_case_folder/drep.csv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --pathways "$use_case_folder/kegg_pathway_completeness.tsv" \
  --top_modules None \
  --mode mean \
  --output out
  ```

Available modes:
- `mean`: highlight core functional modules
- `variance`: highlight modules differing most between clusters
- `both`: create both plots

### Pathway Module Heatmap

Instead of the drep cluster plot above, it is possible to produce just the heatmap without the bars.

```bash
mags-visualization pathway-module-heatmap \
  --drep "$use_case_folder/drep.csv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --pathways "$use_case_folder/kegg_pathway_completeness.tsv" \
  --top_modules 30 \
  --top_representatives 30 \
  --representatives_only \
  --sort_by cluster_size \
  --show_module_labels \
  --output out
```

Use `--show_module_labels` to display KEGG module accessions on the x-axis.
For larger plots, omitting labels can improve readability.

### Completeness / Contamination

```bash
mags-visualization comp-conta \
  --checkm "$use_case_folder/checkm.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --mode quality \
  --output out
```

### Taxonomy Sankey plot

```bash
mags-visualization taxa-sankey \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --rank phylum \
  --output out
```

## Macroalgal microbiome

```bash
use_case_folder="marine-use-case/data"
```

### Sample heatmap

* Family

```bash
mags-visualization sample-heatmap \
  --coverm "$use_case_folder/coverm.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --tax_level family \
  --no_log \
  --top_bar_height 2.0 \
  --top_bar_spacer -1 \
  --output "$use_case_folder/plots"
```

* Species

```bash
mags-visualization sample-heatmap \
  --coverm "$use_case_folder/coverm.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --tax_level species \
  --output "$use_case_folder/plots"
```

### dRep cluster plot

```bash
mags-visualization drep-cluster-annot \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --drep "$use_case_folder/drep.csv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --bakta "$use_case_folder/bakta.tsv" \
  --quast "$use_case_folder/quast.tsv" \
  --output out
```

### dRep cluster plot with functional annotation (KEGG pathway completeness)

```bash
mags-visualization drep-cluster-func \
  --drep "$use_case_folder/drep.csv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --pathways "$use_case_folder/kegg_pathway_completeness.tsv" \
  --top_modules 30 \
  --mode mean \
  --output out
  ```

To include all modules, use `--top_modules None`.

```bash
mags-visualization drep-cluster-func \
  --drep "$use_case_folder/drep.csv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --pathways "$use_case_folder/kegg_pathway_completeness.tsv" \
  --top_modules None \
  --mode mean \
  --output out
  ```

### Pathway Module Heatmap

```bash
mags-visualization pathway-module-heatmap \
  --drep "$use_case_folder/drep.csv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --pathways "$use_case_folder/kegg_pathway_completeness.tsv" \
  --top_modules 30 \
  --top_representatives 30 \
  --representatives_only \
  --sort_by cluster_size \
  --show_module_labels \
  --output out
```

### Completeness / Contamination

```bash
mags-visualization comp-conta \
  --checkm "$use_case_folder/checkm.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --mode quality \
  --output out
```

### Taxonomy Sankey plot

```bash
mags-visualization taxa-sankey \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --rank phylum \
  --output out
```


## Cloud

```bash
use_case_folder="cloud-use-case/data"
```

### Sample heatmap

```bash
mags-visualization sample-heatmap \
  --coverm "$use_case_folder/coverm.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --metadata "$use_case_folder/metadata.tsv" \
  --meta_cols "Geographic origin of the air mass" "Condition" "Season" \
  --tax_level genus \
  --top_bar_spacer -0.5 \
  --no_log \
  --top_bar_height 3.5 \
  --output out
```

If without --no_log:
```bash
--top_bar_height 1.4
```

### dRep cluster plot

```bash
mags-visualization drep-cluster-annot \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --drep "$use_case_folder/drep.csv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --bakta "$use_case_folder/bakta.tsv" \
  --quast "$use_case_folder/quast.tsv" \
  --output out
```

### dRep cluster plot with functional annotation (KEGG pathway completeness)

```bash
mags-visualization drep-cluster-func \
  --drep "$use_case_folder/drep.csv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --pathways "$use_case_folder/kegg_pathway_completeness.tsv" \
  --top_modules 30 \
  --mode mean \
  --output out
  ```

To include all modules, use `--top_modules None`.

```bash
mags-visualization drep-cluster-func \
  --drep "$use_case_folder/drep.csv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --pathways "$use_case_folder/kegg_pathway_completeness.tsv" \
  --top_modules None \
  --mode mean \
  --output out
  ```

### Pathway Module Heatmap

```bash
mags-visualization pathway-module-heatmap \
  --drep "$use_case_folder/drep.csv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --pathways "$use_case_folder/kegg_pathway_completeness.tsv" \
  --top_modules 30 \
  --top_representatives 30 \
  --representatives_only \
  --sort_by cluster_size \
  --show_module_labels \
  --output out
```

### Completeness / Contamination

```bash
mags-visualization comp-conta \
  --checkm "$use_case_folder/checkm.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --mode quality \
  --output out
```

### Taxonomy Sankey plot

```bash
mags-visualization taxa-sankey \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --rank phylum \
  --output out
```

## Termite Head Microbiome

```bash
use_case_folder="termite-use-case/data"
```

### Sample heatmap

* Genus

```bash
mags-visualization sample-heatmap \
  --coverm "$use_case_folder/coverm.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --metadata "$use_case_folder/metadata.tsv" \
  --meta_cols "Species" "Casts" "Colony" \
  --tax_level genus \
  --no_log \
  --spacer_meta 4.5 \
  --top_bar_height 1.5 \
  --hspace 0.08 \
  --output "$use_case_folder/plots"
```

* Species

```bash
mags-visualization sample-heatmap \
  --coverm "$use_case_folder/coverm.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --metadata "$use_case_folder/metadata.tsv" \
  --meta_cols "Species" "Casts" "Colony" \
  --tax_level species \
  --no_log \
  --spacer_meta 4.5 \
  --top_bar_height 1.5 \
  --hspace 0.08 \
  --output "$use_case_folder/plots"
```

* Family

```bash
mags-visualization sample-heatmap \
  --coverm "$use_case_folder/coverm.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --metadata "$use_case_folder/metadata.tsv" \
  --meta_cols "Species" "Casts" "Colony" \
  --tax_level family \
  --no_log \
  --spacer_meta 4.5 \
  --top_bar_height 1.5 \
  --hspace 0.08 \
  --output "$use_case_folder/plots"
```

### dRep cluster plot

```bash
mags-visualization drep-cluster-annot \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --drep "$use_case_folder/drep.csv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --quast "$use_case_folder/quast.tsv" \
  --output out
```

### Completeness / Contamination

```bash
mags-visualization comp-conta \
  --checkm "$use_case_folder/checkm.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --mode quality \
  --output out
```

### Taxonomy Sankey plot

```bash
mags-visualization taxa-sankey \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --rank phylum \
  --output out
```


# Powershell
If you want to use powershell instead of bash, it's just a different syntax.
In the following examples the **subcommand** "all" is used, which produces all plots.

## Macroalgal microbiome

```powershell
$use_case_folder = "marine-use-case\data"

mags-visualization all `
  --coverm "$use_case_folder\coverm.tsv" `
  --checkm "$use_case_folder\checkm.tsv" `
  --checkm2 "$use_case_folder\checkm2.tsv" `
  --gtdb "$use_case_folder\gtdb.tsv" `
  --drep "$use_case_folder\drep.csv" `
  --quast "$use_case_folder\quast.tsv" `
  --bakta "$use_case_folder\bakta.tsv" `
  --pathways "$use_case_folder\kegg_pathway_completeness.tsv" `
  --color_by tax `
  --tax_level genus `
  --max_col 10 `
  --top_n 30 `
  -o test
```

## Bee

```powershell
$use_case_folder = "bee-use-case\data"

mags-visualization all `
  --coverm "$use_case_folder\coverm.tsv" `
  --checkm "$use_case_folder\checkm.tsv" `
  --checkm2 "$use_case_folder\checkm2.tsv" `
  --gtdb "$use_case_folder\gtdb.tsv" `
  --drep "$use_case_folder\drep.csv" `
  --quast "$use_case_folder\quast.tsv" `
  --bakta "$use_case_folder\bakta.tsv" `
  --pathways "$use_case_folder\kegg_pathway_completeness.tsv" `
  --metadata "$use_case_folder\metadata.tsv" `
  --meta_cols "Infection by Nosema ceranae" "Chronic exposure to neonicotinoid" "Treatment with probiotic" `
  --color_by tax `
  --tax_level phylum `
  --top_n 30 `
  --top_bar_spacer -0.5 `
  --spacer_meta 2.5 `
  -o test
```

```powershell
--no_log `
--top_bar_height 2.0
```

## Cloud

```powershell
$use_case_folder = "cloud-use-case\data"

mags-visualization all `
  --coverm "$use_case_folder\coverm.tsv" `
  --checkm "$use_case_folder\checkm.tsv" `
  --checkm2 "$use_case_folder\checkm2.tsv" `
  --gtdb "$use_case_folder\gtdb.tsv" `
  --drep "$use_case_folder\drep.csv" `
  --quast "$use_case_folder\quast.tsv" `
  --bakta "$use_case_folder\bakta.tsv" `
  --pathways "$use_case_folder\kegg_pathway_completeness.tsv" `
  --metadata "$use_case_folder\metadata.tsv" `
  --meta_cols "Geographic origin of the air mass" "Condition" "Season" `
  --color_by tax `
  --tax_level phylum `
  --top_n 30 `
  --top_bar_spacer -0.5 `
  --no_log `
  --top_bar_height 3.5 `
  -o test
```

If without --no_log:
```powershell
--top_bar_height 1.4
```

## Termite

```powershell
$use_case_folder = "termite-use-case\data"

mags-visualization all `
  --coverm "$use_case_folder\coverm.tsv" `
  --checkm "$use_case_folder\checkm.tsv" `
  --checkm2 "$use_case_folder\checkm2.tsv" `
  --gtdb "$use_case_folder\gtdb.tsv" `
  --drep "$use_case_folder\drep.csv" `
  --quast "$use_case_folder\quast.tsv" `
  --metadata "$use_case_folder\metadata.tsv" `
  --meta_cols "Species" "Casts" "Colony" `
  --color_by tax `
  --tax_level genus `
  --max_col 10 `
  --top_n 30 `
  --spacer_meta 4.0 `
  --no_log `
  -o test
```
