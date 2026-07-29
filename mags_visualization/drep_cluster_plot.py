import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec, cm
from matplotlib.colors import Normalize
from .id_normalizer import normalize_genome_id, format_genome_display

# ---- Tax and ID processing ----
def extract_rank(tax, rank: str):
    """Extract a GTDB rank from a taxonomy string."""
    if pd.isna(tax):
        return "Unclassified"
    
    prefix_map = {
        "domain": "d__",
        "phylum": "p__",
        "class": "c__",
        "order": "o__",
        "family": "f__",
        "genus": "g__",
        "species": "s__",
    }
    
    prefix = prefix_map.get(rank.lower())
    if not prefix:
        raise ValueError(f"Invalid rank '{rank}'. Choose one of {list(prefix_map.keys())}.")
    
    for part in str(tax).split(";"):
        part = part.strip()
        if part.startswith(prefix):
            name = part[len(prefix):].strip()
            return name if name else "Unclassified"
    
    return "Unclassified"


# ---- Main plot function ----
def drep_cluster_plot(drep_df: pd.DataFrame, gtdb_df: pd.DataFrame, output_path: str, 
                     tax_levels=("phylum", "genus"), top_n: int = 30, fmt: str = "png",
                     tax_levels_space: float = 0.3, fig_size=None, checkm2_df: pd.DataFrame = None,
                     quast_df: pd.DataFrame = None, bakta_df: pd.DataFrame = None,
                     require_quality: bool = False,):
    """
    Create a horizontal bar plot showing the top N clusters by member count,
    annotated with:
    - number of MAGs
    - representative ID
    - representative taxonomy - phylum and genus
    - heatmap:
        - completeness and contamination
        - genome size and GC % from QUAST
        - CDS, rRNA counts and hypotheticals/CDS ratio
    """
    
    os.makedirs(output_path, exist_ok=True)

    # ---- Validate tax levels ----
    valid_ranks = {"domain", "phylum", "class", "order", "family", "genus", "species"}
    if tax_levels is None:
        tax_levels = ("phylum", "genus")

    if isinstance(tax_levels, str):
        tax_levels = tuple(tax_levels.split())

    tax_levels = tuple([t.lower().strip() for t in tax_levels if str(t).strip() != ""])
    if len(tax_levels) == 0:
        tax_levels = ("phylum", "genus")

    for lvl in tax_levels:
        if lvl not in valid_ranks:
            raise ValueError(f"Invalid tax level '{lvl}'. Choose from {sorted(valid_ranks)}.")

    print(f"[INFO] Creating dRep cluster plot (top {top_n}, tax_levels={tax_levels})")

    # ---- dRep & gtdb ----
    drep = drep_df.copy()
    gtdb = gtdb_df.copy()
    
    
    if drep.index.name is not None:
        drep = drep.reset_index()
    
    # Normalize genome IDs
    if 'genome' in drep.columns:
        drep['genome_normalized'] = drep['genome'].apply(normalize_genome_id)
    else:
        raise ValueError("dRep DataFrame must contain 'genome' column")
    
    if 'secondary_cluster' not in drep.columns:
        raise ValueError("dRep DataFrame must contain 'secondary_cluster' column")
    
    # Process GTDB data
    if gtdb.index.name is not None:
        gtdb = gtdb.reset_index()

    if 'user_genome' not in gtdb.columns or 'classification' not in gtdb.columns:
        raise ValueError("GTDB file must contain 'user_genome' and 'classification' columns")

    gtdb['user_genome_normalized'] = gtdb['user_genome'].apply(normalize_genome_id)

    # Extract requested ranks
    for lvl in set(tax_levels):
        gtdb[lvl] = gtdb["classification"].apply(lambda x: extract_rank(x, lvl))

    gtdb_genomes = set(gtdb['user_genome_normalized'].values)

    def get_representative(cluster_genomes):
        """Find the genome from this cluster that is in GTDB (the representative)"""
        for genome in cluster_genomes:
            if genome in gtdb_genomes:
                return genome
        # if no genome found in GTDB, take first one
        return cluster_genomes.iloc[0]

    cluster_reps = drep.groupby('secondary_cluster')['genome_normalized'].apply(get_representative).reset_index()
    cluster_reps.columns = ['secondary_cluster', 'representative']

    print(f"[INFO] Found {len(cluster_reps)} representatives, {cluster_reps['representative'].isin(gtdb_genomes).sum()} are in GTDB")

    # ---- Optional: filter to clusters with CheckM2 quality data ----
    if require_quality:
        if checkm2_df is None:
            print("[WARN] require_quality=True but no checkm2_df provided; skipping filter.")
        else:
            checkm2 = checkm2_df.copy()
            if checkm2.index.name is not None:
                checkm2 = checkm2.reset_index()
                checkm2.rename(columns={checkm2.columns[0]: 'Name'}, inplace=True)
            if 'Name' in checkm2.columns:
                checkm2['genome_normalized'] = checkm2['Name'].apply(normalize_genome_id)
                checkm2_genomes = set(checkm2['genome_normalized'].values)
                before = len(cluster_reps)
                cluster_reps = cluster_reps[cluster_reps['representative'].isin(checkm2_genomes)]
                print(f"[INFO] require_quality: {before} -> {len(cluster_reps)} clusters after CheckM2 filter")
            else:
                print("[WARN] CheckM2 DataFrame has no 'Name' column; cannot filter by quality.")
        if len(cluster_reps) == 0:
            print("[WARN] No clusters pass the quality filter; returning empty plot.")

    # Count members per cluster (after optional filter)
    cluster_counts = drep.groupby('secondary_cluster').size().reset_index(name='n_members')
    cluster_counts = cluster_counts[cluster_counts['secondary_cluster'].isin(cluster_reps['secondary_cluster'])]
    cluster_counts = cluster_counts.sort_values('n_members', ascending=False).head(top_n)

    total_clusters = drep['secondary_cluster'].nunique()
    filtered_clusters = len(cluster_reps)
    if require_quality:
        print(f"[INFO] Found {total_clusters} total secondary clusters, filtered to {filtered_clusters} with QC data, showing top {top_n}")
    else:
        print(f"[INFO] Found {total_clusters} total secondary clusters, showing top {top_n}")

    cluster_data = cluster_counts.merge(cluster_reps, on='secondary_cluster', how='left')
    cluster_data['representative_display'] = cluster_data['representative'].apply(format_genome_display)
    
    cols_to_merge = ['user_genome_normalized'] + list(tax_levels)

    # Attach tax to cluster data
    cluster_data = cluster_data.merge(
        gtdb[cols_to_merge],
        left_on='representative',
        right_on='user_genome_normalized',
        how='left'
    )

    # Fill taxonomy
    for lvl in tax_levels:
        if lvl in cluster_data.columns:
            cluster_data[lvl] = cluster_data[lvl].fillna("Unclassified")
        else:
            cluster_data[lvl] = "Unclassified"

    cluster_data = cluster_data.sort_values("n_members", ascending=True)

    # ---- CheckM2 data - completeness and contamination ----
    rep_quality_available = False
    if checkm2_df is not None:
        print("[INFO] Processing CheckM2 data for representative completeness/contamination")
        checkm2 = checkm2_df.copy()
        
        if checkm2.index.name is not None:
            checkm2 = checkm2.reset_index()
            checkm2.rename(columns={checkm2.columns[0]: 'Name'}, inplace=True)
        
        # Normalize genome IDs in CheckM2
        if 'Name' in checkm2.columns:
            checkm2['genome_normalized'] = checkm2['Name'].apply(normalize_genome_id)
        else:
            print("[WARN] CheckM2 DataFrame has no 'Name' column, skipping quality")
            checkm2 = None
        
        if checkm2 is not None:
            cluster_data = cluster_data.merge(
                checkm2[['genome_normalized', 'Completeness', 'Contamination']],
                left_on='representative',
                right_on='genome_normalized',
                how='left',
                suffixes=('', '_checkm2')
            )
            cluster_data = cluster_data.drop(columns=['genome_normalized'], errors='ignore')
            cluster_data = cluster_data.rename(columns={
                'Completeness': 'rep_completeness',
                'Contamination': 'rep_contamination'
            })
            rep_quality_available = True
        
    # ---- QUAST -> representative genome_size, gc ----
    quast_available = False
    if quast_df is not None:
        quast = quast_df.copy()

        gc_row_name = None
        size_row_name = None

        for name in quast.index:
            lname = str(name).lower()
            if gc_row_name is None and 'gc' in lname and '%' in lname:
                gc_row_name = name
            if size_row_name is None and 'total length' in lname:
                size_row_name = name

        if gc_row_name is not None and size_row_name is not None:
            gc_series = quast.loc[gc_row_name]
            size_series = quast.loc[size_row_name]

            quast_long = pd.DataFrame({
                'genome_raw': gc_series.index,
                'gc': pd.to_numeric(gc_series.values, errors='coerce'),
                'genome_size': pd.to_numeric(size_series.values, errors='coerce'),
            })

            quast_long['genome_normalized'] = quast_long['genome_raw'].apply(normalize_genome_id)

            # Merge with cluster_data
            cluster_data = cluster_data.merge(
                quast_long[['genome_normalized', 'gc', 'genome_size']],
                left_on='representative',
                right_on='genome_normalized',
                how='left',
                suffixes=('', '_quast')
            )
            
            cluster_data = cluster_data.rename(columns={
                'gc': 'rep_gc',
                'genome_size': 'rep_genome_size'
            })
            quast_available = True
        else:
            print("[WARN] Could not find GC (%) or Total length rows in QUAST file; skipping GC/genome size.")
    
    # ---- Bakta -> representative cds, rrna, hypo_cds_ratio ----
    bakta_available = False
    if bakta_df is not None:
        bakta = bakta_df.copy()

        def _find_row(name_substr):
            name_substr = name_substr.lower()
            for idx in bakta.index:
                if name_substr in str(idx).lower():
                    return idx
            return None

        cds_row_name  = _find_row("cds")
        rrna_row_name = _find_row("rrna")
        hypo_row_name = _find_row("hypothetical")

        if cds_row_name is not None and rrna_row_name is not None:
            cds_series  = bakta.loc[cds_row_name]
            rrna_series = bakta.loc[rrna_row_name]
            hypo_series = bakta.loc[hypo_row_name] if hypo_row_name is not None else None

            bakta_long = pd.DataFrame({
                "genome_raw": cds_series.index,
                "cds":  pd.to_numeric(cds_series.values, errors="coerce"),
                "rrna": pd.to_numeric(rrna_series.values, errors="coerce"),
            })
            
            has_hypotheticals = False
            if hypo_series is not None:
                bakta_long["hypotheticals"] = pd.to_numeric(hypo_series.values, errors="coerce")
                # Calculate ratio at MAG level
                bakta_long["hypo_cds_ratio"] = bakta_long["hypotheticals"] / bakta_long["cds"]
                bakta_long["hypo_cds_ratio"] = bakta_long["hypo_cds_ratio"].replace([np.inf, -np.inf], np.nan).fillna(0)
                has_hypotheticals = True

            bakta_long["genome_raw_clean"] = bakta_long["genome_raw"].str.replace("_Count", "", regex=False)
            bakta_long["genome_normalized"] = bakta_long["genome_raw_clean"].apply(normalize_genome_id)

            merge_cols = ["genome_normalized", "cds", "rrna"]
            if has_hypotheticals:
                merge_cols.extend(["hypotheticals", "hypo_cds_ratio"])
            
            cluster_data = cluster_data.merge(
                bakta_long[merge_cols],
                left_on='representative',
                right_on='genome_normalized',
                how='left',
                suffixes=('', '_bakta')
            )
            
            rename_dict = {"cds": "rep_cds", "rrna": "rep_rrna"}
            if has_hypotheticals:
                rename_dict["hypotheticals"] = "rep_hypotheticals"
                rename_dict["hypo_cds_ratio"] = "rep_hypo_cds_ratio"
            
            cluster_data = cluster_data.rename(columns=rename_dict)
            bakta_available = True
        else:
            print("[WARN] Could not find CDS or rRNA rows in Bakta file; skipping Bakta columns.")

    # ---- Figure & Layout ----
    if fig_size is None:
        if rep_quality_available:
            fig_size = (14, max(7, top_n * 0.28))
        else:
            fig_size = (14, max(8, top_n * 0.3))

    fig = plt.figure(figsize=fig_size)

    n_tax = len(tax_levels)

    if rep_quality_available:
        # [sample names | bars | spacer | tax... | spacer | heatmap]
        width_ratios = [0.03, 2, 0.02] + [tax_levels_space] * n_tax + [0.04, 0.4]
        gs = gridspec.GridSpec(
            1,
            len(width_ratios),
            figure=fig,
            width_ratios=width_ratios,
            wspace=0.04,
        )
        ax_samples = fig.add_subplot(gs[0, 0])
        ax_bars = fig.add_subplot(gs[0, 1])

        tax_axes = []
        tax_start_col = 3
        for i, lvl in enumerate(tax_levels):
            tax_axes.append((lvl, fig.add_subplot(gs[0, tax_start_col + i])))

        ax_heatmap = fig.add_subplot(gs[0, len(width_ratios) - 1])
    else:
        # [sample names | bars | spacer | tax...]
        width_ratios = [0.15, 2, 0.05] + [0.30] * n_tax
        gs = gridspec.GridSpec(
            1,
            len(width_ratios),
            figure=fig,
            width_ratios=width_ratios,
            wspace=0.08,
        )
        ax_samples = fig.add_subplot(gs[0, 0])
        ax_bars = fig.add_subplot(gs[0, 1])

        tax_axes = []
        tax_start_col = 3
        for i, lvl in enumerate(tax_levels):
            tax_axes.append((lvl, fig.add_subplot(gs[0, tax_start_col + i])))

        ax_heatmap = None

    # ---- Bar Plot ----
    y_pos = np.arange(len(cluster_data))
    
    bars = ax_bars.barh(
        y_pos,
        cluster_data['n_members'].values,
        height=0.7,
        color='#1f2937',
        edgecolor='#111827',
        linewidth=0.5
    )
    
    max_members = int(cluster_data['n_members'].max())
    nice_max = int(np.ceil(max_members * 1.05))
    step = 5 if nice_max <= 50 else 10

    ax_bars.set_xlim(0, nice_max)
    ax_bars.set_xticks(np.arange(0, nice_max + 1, step))
    ax_bars.set_ylim(-0.5, len(cluster_data) - 0.5)
    ax_bars.set_xlabel('# of MAGs in species-level cluster', fontsize=11)
    
    ax_bars.grid(axis='x', color='#e5e7eb', linewidth=0.5, zorder=0)
    ax_bars.set_axisbelow(True)
    
    ax_bars.set_yticks([])
    ax_bars.spines['left'].set_visible(False)
    ax_bars.spines['top'].set_visible(False)
    ax_bars.spines['right'].set_visible(False)
    
    for i, (idx, row) in enumerate(cluster_data.iterrows()):
        count = int(row['n_members'])
        ax_bars.text(
            count + nice_max * 0.01,
            i,
            str(count),
            va='center',
            ha='left',
            fontsize=9,
            fontweight='bold'
        )
    
    # ---- Configure sample names axis ----
    ax_samples.set_xlim(0, 1)
    ax_samples.set_ylim(-0.5, len(cluster_data) - 0.5)
    ax_samples.set_yticks(y_pos)
    ax_samples.set_yticklabels(cluster_data['representative_display'].values, fontsize=8, ha='right')
    ax_samples.tick_params(
        axis="y",
        which="both",
        length=10,
        width=0.4,
        pad=4
    )
    ax_samples.set_xticks([])

    for spine in ax_samples.spines.values():
        spine.set_visible(False)

    # ---- Configure taxonomy axes (dynamic) ----
    for lvl, ax in tax_axes:
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.5, len(cluster_data) - 0.5)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(cluster_data[lvl].values, fontsize=8, ha="left")
        ax.set_xticks([])
        ax.tick_params(axis="y", which="both", length=0, pad=2)

        for spine in ax.spines.values():
            spine.set_visible(False)

        ax.text(
            0.0,
            1.02,
            lvl.capitalize(),
            ha="left",
            va="bottom",
            fontsize=9,
            fontweight="bold",
            transform=ax.transAxes,
        )

        ax.grid(False)
        ax.xaxis.grid(False)
        ax.yaxis.grid(False)
        
    # ---- Heatmap - representative values only ----
    if rep_quality_available and ax_heatmap is not None:

        ax_heatmap.set_ylim(-0.5, len(cluster_data) - 0.5)
        ax_heatmap.set_yticks([])
        ax_heatmap.spines['left'].set_visible(False)
        ax_heatmap.spines['top'].set_visible(False)
        ax_heatmap.spines['right'].set_visible(False)

        block_width = 1.0
        block_pos = {}
        x = 0

        # Completeness and Contamination blocks
        block_pos['completeness'] = x
        block_pos['contamination'] = x + 1
        x += 2

        if quast_available:
            block_pos['size'] = x
            block_pos['gc']   = x + 1
            x += 2

        if bakta_available:
            block_pos['cds']  = x
            block_pos['rrna'] = x + 1
            x += 2
            if 'rep_hypo_cds_ratio' in cluster_data.columns:
                block_pos['hypo_ratio'] = x
                x += 1

        n_blocks = x
        ax_heatmap.set_xlim(0, n_blocks)
        ax_heatmap.set_xticks([])

        # --- Normalizing + Colormaps ---
        comp_norm = cont_norm = size_norm = gc_norm = cds_norm = rrna_norm = hypo_ratio_norm = None
        comp_cmap = cont_cmap = size_cmap = gc_cmap = cds_cmap = rrna_cmap = hypo_ratio_cmap = None

        # Completeness and Contamination
        if rep_quality_available:
            comp_vals = cluster_data['rep_completeness'].values.astype(float)
            cont_vals = cluster_data['rep_contamination'].values.astype(float)
            comp_vals_valid = comp_vals[~np.isnan(comp_vals)]
            cont_vals_valid = cont_vals[~np.isnan(cont_vals)]
            
            if len(comp_vals_valid) > 0:
                comp_norm = Normalize(vmin=comp_vals_valid.min(), vmax=comp_vals_valid.max())
                comp_cmap = plt.get_cmap('Greens')
            
            if len(cont_vals_valid) > 0:
                cont_norm = Normalize(vmin=cont_vals_valid.min(), vmax=cont_vals_valid.max())
                cont_cmap = plt.get_cmap('Reds')

        # Quast (genome size, GC%)
        if quast_available:
            size_vals = cluster_data['rep_genome_size'].values.astype(float)
            gc_vals   = cluster_data['rep_gc'].values.astype(float)
            size_vals_valid = size_vals[~np.isnan(size_vals)]
            gc_vals_valid   = gc_vals[~np.isnan(gc_vals)]
            if len(size_vals_valid) > 0:
                size_norm = Normalize(vmin=size_vals_valid.min(), vmax=size_vals_valid.max())
                size_cmap = plt.get_cmap('Blues')
            if len(gc_vals_valid) > 0:
                gc_norm = Normalize(vmin=gc_vals_valid.min(), vmax=gc_vals_valid.max())
                gc_cmap = plt.get_cmap('Purples')

        # Bakta (CDS, rRNA, hypotheticals ratio)
        if bakta_available:
            cds_vals  = cluster_data['rep_cds'].values.astype(float)
            rrna_vals = cluster_data['rep_rrna'].values.astype(float)
            cds_vals_valid  = cds_vals[~np.isnan(cds_vals)]
            rrna_vals_valid = rrna_vals[~np.isnan(rrna_vals)]
            if len(cds_vals_valid) > 0:
                cds_norm = Normalize(vmin=cds_vals_valid.min(), vmax=cds_vals_valid.max())
                cds_cmap = plt.get_cmap('Oranges')
            if len(rrna_vals_valid) > 0:
                rrna_norm = Normalize(vmin=rrna_vals_valid.min(), vmax=rrna_vals_valid.max())
                rrna_cmap = plt.get_cmap('YlOrBr')
            
            if 'rep_hypo_cds_ratio' in cluster_data.columns:
                hypo_ratio_vals = cluster_data['rep_hypo_cds_ratio'].values.astype(float)
                hypo_ratio_vals_valid = hypo_ratio_vals[~np.isnan(hypo_ratio_vals)]
                if len(hypo_ratio_vals_valid) > 0:
                    hypo_ratio_norm = Normalize(vmin=hypo_ratio_vals_valid.min(), vmax=hypo_ratio_vals_valid.max())
                    hypo_ratio_cmap = plt.get_cmap('RdPu')

        # Draw a row of colored blocks per cluster
        for i, row in cluster_data.iterrows():
            # Completeness
            if rep_quality_available and comp_norm is not None:
                comp = row['rep_completeness']
                col_comp = comp_cmap(comp_norm(comp)) if pd.notna(comp) else 'white'
                ax_heatmap.barh(
                    i, block_width, left=block_pos['completeness'],
                    height=0.7, color=col_comp, edgecolor='white', linewidth=0.2
                )
            
            # Contamination
            if rep_quality_available and cont_norm is not None:
                cont = row['rep_contamination']
                col_cont = cont_cmap(cont_norm(cont)) if pd.notna(cont) else 'white'
                ax_heatmap.barh(
                    i, block_width, left=block_pos['contamination'],
                    height=0.7, color=col_cont, edgecolor='white', linewidth=0.2
                )

            # QUAST
            if quast_available and size_norm is not None:
                size = row['rep_genome_size']
                col_size = size_cmap(size_norm(size)) if pd.notna(size) else 'white'
                ax_heatmap.barh(
                    i, block_width, left=block_pos['size'],
                    height=0.7, color=col_size, edgecolor='white', linewidth=0.2
                )

            if quast_available and gc_norm is not None:
                gc_val = row['rep_gc']
                col_gc = gc_cmap(gc_norm(gc_val)) if pd.notna(gc_val) else 'white'
                ax_heatmap.barh(
                    i, block_width, left=block_pos['gc'],
                    height=0.7, color=col_gc, edgecolor='white', linewidth=0.2
                )

            # Bakta
            if bakta_available and cds_norm is not None:
                cds = row['rep_cds']
                col_cds = cds_cmap(cds_norm(cds)) if pd.notna(cds) else 'white'
                ax_heatmap.barh(
                    i, block_width, left=block_pos['cds'],
                    height=0.7, color=col_cds, edgecolor='white', linewidth=0.2
                )

            if bakta_available and rrna_norm is not None:
                rr = row['rep_rrna']
                col_rr = rrna_cmap(rrna_norm(rr)) if pd.notna(rr) else 'white'
                ax_heatmap.barh(
                    i, block_width, left=block_pos['rrna'],
                    height=0.7, color=col_rr, edgecolor='white', linewidth=0.2
                )
            
            if bakta_available and 'hypo_ratio' in block_pos and hypo_ratio_norm is not None:
                ratio = row['rep_hypo_cds_ratio']
                col_ratio = hypo_ratio_cmap(hypo_ratio_norm(ratio)) if pd.notna(ratio) else 'white'
                ax_heatmap.barh(
                    i, block_width, left=block_pos['hypo_ratio'],
                    height=0.7, color=col_ratio, edgecolor='white', linewidth=0.2
                )

        # Vertical labels over each block
        label_y = len(cluster_data) - 0.5 + 0.5
        
        ax_heatmap.text(
            block_pos['completeness'] + 0.5, label_y, 'Compl.',
            ha='center', va='bottom', fontsize=8, rotation=90
        )
        ax_heatmap.text(
            block_pos['contamination'] + 0.5, label_y, 'Cont.',
            ha='center', va='bottom', fontsize=8, rotation=90
        )
        
        if quast_available:
            ax_heatmap.text(
                block_pos['size'] + 0.5, label_y, 'Size',
                ha='center', va='bottom', fontsize=8, rotation=90
            )
            ax_heatmap.text(
                block_pos['gc'] + 0.5, label_y, 'GC',
                ha='center', va='bottom', fontsize=8, rotation=90
            )
        
        if bakta_available:
            ax_heatmap.text(
                block_pos['cds'] + 0.5, label_y, 'CDS',
                ha='center', va='bottom', fontsize=8, rotation=90
            )
            ax_heatmap.text(
                block_pos['rrna'] + 0.5, label_y, 'rRNA',
                ha='center', va='bottom', fontsize=8, rotation=90
            )
            
            if 'hypo_ratio' in block_pos:
                ax_heatmap.text(
                    block_pos['hypo_ratio'] + 0.5, label_y, 'Hypo/\nCDS',
                    ha='center', va='bottom', fontsize=8, rotation=90
                )

        pos = ax_heatmap.get_position()
        legend_y = pos.y0 - 0.15

        cbar_h = 0.10
        cbar_w = 0.015
        spacing = 0.04

        current_x = pos.x0 - 0.25
        
        # Colorbars
        if rep_quality_available and comp_norm is not None:
            ax_cbar_comp = fig.add_axes([current_x, legend_y, cbar_w, cbar_h])
            comp_cbar = cm.ScalarMappable(norm=comp_norm, cmap=comp_cmap)
            comp_cbar.set_array([])
            cb1 = fig.colorbar(comp_cbar, cax=ax_cbar_comp, orientation="vertical")
            cb1.ax.tick_params(labelsize=6)
            cb1.set_label("Completeness (%)", fontsize=7, labelpad=1)
            current_x += cbar_w + spacing

        if rep_quality_available and cont_norm is not None:
            ax_cbar_cont = fig.add_axes([current_x, legend_y, cbar_w, cbar_h])
            cont_cbar = cm.ScalarMappable(norm=cont_norm, cmap=cont_cmap)
            cont_cbar.set_array([])
            cb2 = fig.colorbar(cont_cbar, cax=ax_cbar_cont, orientation="vertical")
            cb2.ax.tick_params(labelsize=6)
            cb2.set_label("Contamination (%)", fontsize=7, labelpad=1)
            current_x += cbar_w + spacing

        if quast_available and size_norm is not None:
            ax_cbar_size = fig.add_axes([current_x, legend_y, cbar_w, cbar_h])
            size_cbar = cm.ScalarMappable(norm=size_norm, cmap=size_cmap)
            size_cbar.set_array([])
            cb3 = fig.colorbar(size_cbar, cax=ax_cbar_size, orientation="vertical")
            cb3.ax.tick_params(labelsize=6)
            cb3.set_label("Genome size (bp)", fontsize=7, labelpad=1)
            current_x += cbar_w + spacing

        if quast_available and gc_norm is not None:
            ax_cbar_gc = fig.add_axes([current_x, legend_y, cbar_w, cbar_h])
            gc_cbar = cm.ScalarMappable(norm=gc_norm, cmap=gc_cmap)
            gc_cbar.set_array([])
            cb4 = fig.colorbar(gc_cbar, cax=ax_cbar_gc, orientation="vertical")
            cb4.ax.tick_params(labelsize=6)
            cb4.set_label("GC (%)", fontsize=7, labelpad=1)
            current_x += cbar_w + spacing

        if bakta_available and cds_norm is not None:
            ax_cbar_cds = fig.add_axes([current_x, legend_y, cbar_w, cbar_h])
            cds_cbar = cm.ScalarMappable(norm=cds_norm, cmap=cds_cmap)
            cds_cbar.set_array([])
            cb5 = fig.colorbar(cds_cbar, cax=ax_cbar_cds, orientation="vertical")
            cb5.ax.tick_params(labelsize=6)
            cb5.set_label("CDS", fontsize=7, labelpad=1)
            current_x += cbar_w + spacing

        if bakta_available and rrna_norm is not None:
            ax_cbar_rrna = fig.add_axes([current_x, legend_y, cbar_w, cbar_h])
            rrna_cbar = cm.ScalarMappable(norm=rrna_norm, cmap=rrna_cmap)
            rrna_cbar.set_array([])
            cb6 = fig.colorbar(rrna_cbar, cax=ax_cbar_rrna, orientation="vertical")
            cb6.ax.tick_params(labelsize=6)
            cb6.set_label("rRNAs", fontsize=7, labelpad=1)
            current_x += cbar_w + spacing
        
        if bakta_available and hypo_ratio_norm is not None:
            ax_cbar_hypo = fig.add_axes([current_x, legend_y, cbar_w, cbar_h])
            hypo_cbar = cm.ScalarMappable(norm=hypo_ratio_norm, cmap=hypo_ratio_cmap)
            hypo_cbar.set_array([])
            cb7 = fig.colorbar(hypo_cbar, cax=ax_cbar_hypo, orientation="vertical")
            cb7.ax.tick_params(labelsize=6)
            cb7.set_label("Hypo/CDS ratio", fontsize=7, labelpad=1)

    # ---- Info box ----
    mags_shown = drep[drep['secondary_cluster'].isin(cluster_data['secondary_cluster'])].shape[0]
    clusters_shown = cluster_data['secondary_cluster'].nunique()
    
    info_text = (
        f"Total MAGs: {mags_shown}\n"
        f"Total clusters: {clusters_shown}\n"
        f"MAGs in top {top_n}: {mags_shown} "
        f"({mags_shown/len(drep)*100:.1f}%)"
    )
    
    fig.text(
        0.02, -0.01,
        info_text,
        fontsize=10,
        verticalalignment='bottom',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3, pad=0.5),
        family='monospace'
    )

    plt.suptitle(
        f'Top {top_n} species-level clusters by MAG count',
        fontsize=13,
        y=0.98,
        fontweight='bold'
    )

    txt_file = os.path.join(output_path, "cluster_data_debug.txt")

    with open(txt_file, "w", encoding="utf-8") as f:
        f.write(cluster_data.to_string())

    print(f"[INFO] Debug table written to {txt_file}")
        
    tax_tag = "-".join(tax_levels)
    out_file = os.path.join(output_path, f"drep_cluster_top{top_n}_{tax_tag}.{fmt}")
    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"[INFO] dRep cluster plot saved to {out_file}")
    return out_file