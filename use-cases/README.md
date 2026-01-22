# Example use case
All commands below assume that the command-line tool `mags-visualization`is installed and available in your environment.

The working directory should be **MAGs-visualization/use-cases** 
so that the example paths resolve correctly.

If the command is not available, install the package first (['README.md'](../README.md))

# Commands to generate the use case figure output

## Marine

```bash
use_case_folder="marine-use-case/data"
use_case_output_folder="bee-use-case/plots"

mags-visualization \
  --coverm "$use_case_folder/coverm.tsv" \
  --checkm "$use_case_folder/checkm.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --drep "$use_case_folder/drep.csv" \
  --quast "$use_case_folder/quast.tsv" \
  --bakta "$use_case_folder/bakta.tsv" \
  --color_by tax \
  --tax_level genus \
  --max_col 10 \
  --top_n 30 \
  -o "$use_case_output_folder" \
```

## Bee

### Phylum tax level

```bash
use_case_folder="bee-use-case/data"
use_case_output_folder="bee-use-case/plots"

mags-visualization \
  --coverm "$use_case_folder/coverm.tsv" \
  --checkm "$use_case_folder/checkm.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --drep "$use_case_folder/drep.csv" \
  --quast "$use_case_folder/quast.tsv" \
  --bakta "$use_case_folder/bakta.tsv" \
  --metadata "$use_case_folder/metadata.tsv" \
  --meta_cols "Infection by Nosema ceranae" "Chronic exposure to neonicotinoid" "Treatment with probiotic" \
  --color_by tax \
  --tax_level phylum \
  --top_n 30 \
  --top_bar_spacer -1.5 \
  --spacer_meta 2.5 \
  -o "$use_case_output_folder" \
  --no_log \
  --top_bar_height 3.0
```

### Genus tax level

```bash
use_case_folder="bee-use-case/data"
use_case_output_folder="bee-use-case/plots"

mags-visualization \
  --coverm "$use_case_folder/coverm.tsv" \
  --checkm "$use_case_folder/checkm.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --drep "$use_case_folder/drep.csv" \
  --quast "$use_case_folder/quast.tsv" \
  --bakta "$use_case_folder/bakta.tsv" \
  --metadata "$use_case_folder/metadata.tsv" \
  --meta_cols "Infection by Nosema ceranae" "Chronic exposure to neonicotinoid" "Treatment with probiotic" \
  --color_by tax \
  --tax_level genus \
  --top_n 30 \
  --top_bar_spacer -1.5 \
  --spacer_meta 2.5 \
  -o "$use_case_output_folder" \
  --no_log \
  --top_bar_height 3.0
```

## Cloud

```bash
use_case_folder="cloud-use-case/data"
use_case_output_folder="cloud-use-case/plots"

mags-visualization \
  --coverm "$use_case_folder/coverm.tsv" \
  --checkm "$use_case_folder/checkm.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --drep "$use_case_folder/drep.csv" \
  --quast "$use_case_folder/quast.tsv" \
  --bakta "$use_case_folder/bakta.tsv" \
  --metadata "$use_case_folder/metadata.tsv" \
  --meta_cols "Geographic origin of the air mass" "Condition" "Season" \
  --color_by tax \
  --tax_level phylum \
  --top_n 30 \
  --top_bar_spacer -0.5 \
  --no_log \
  --top_bar_height 3.5 \
  -o "$use_case_output_folder" \
```

If without --no_log:
```bash
--top_bar_height 1.4
```

## Termite

```bash
use_case_folder="termite-use-case/data"
use_case_output_folder="termite-use-case/plots"

mags-visualization \
  --coverm "$use_case_folder/coverm.tsv" \
  --checkm "$use_case_folder/checkm.tsv" \
  --checkm2 "$use_case_folder/checkm2.tsv" \
  --gtdb "$use_case_folder/gtdb.tsv" \
  --drep "$use_case_folder/drep.csv" \
  --quast "$use_case_folder/quast.tsv" \
  --metadata "$use_case_folder/metadata.tsv" \
  --meta_cols "Species" "Casts" "Colony" \
  --color_by tax \
  --tax_level genus \
  --max_col 10 \
  --top_n 30 \
  --spacer_meta 4.0 \
  --no_log \
  -o "$use_case_output_folder" \
```

# Powershell
If you want to use powershell instead of bash, here are the examples.

## Marine

```powershell
$use_case_folder = "marine-use-case\data"
use_case_output_folder="marine-use-case\plots"

mags-visualization `
  --coverm "$use_case_folder\coverm.tsv" `
  --checkm "$use_case_folder\checkm.tsv" `
  --checkm2 "$use_case_folder\checkm2.tsv" `
  --gtdb "$use_case_folder\gtdb.tsv" `
  --drep "$use_case_folder\drep.csv" `
  --quast "$use_case_folder\quast.tsv" `
  --bakta "$use_case_folder\bakta.tsv" `
  --color_by tax `
  --tax_level genus `
  --max_col 10 `
  --top_n 30 `
  -o test
```

## Bee

```powershell
$use_case_folder = "bee-use-case\data"

mags-visualization `
  --coverm "$use_case_folder\coverm.tsv" `
  --checkm "$use_case_folder\checkm.tsv" `
  --checkm2 "$use_case_folder\checkm2.tsv" `
  --gtdb "$use_case_folder\gtdb.tsv" `
  --drep "$use_case_folder\drep.csv" `
  --quast "$use_case_folder\quast.tsv" `
  --bakta "$use_case_folder\bakta.tsv" `
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

mags-visualization `
  --coverm "$use_case_folder\coverm.tsv" `
  --checkm "$use_case_folder\checkm.tsv" `
  --checkm2 "$use_case_folder\checkm2.tsv" `
  --gtdb "$use_case_folder\gtdb.tsv" `
  --drep "$use_case_folder\drep.csv" `
  --quast "$use_case_folder\quast.tsv" `
  --bakta "$use_case_folder\bakta.tsv" `
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

mags-visualization `
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
