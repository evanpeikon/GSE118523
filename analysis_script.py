# analysis_script.py
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from gprofiler.gprofiler import GProfiler

def download_and_prepare_data():
    """Download and prepare data for analysis."""
    os.system("wget -O GSE118523_20161109_old_wt_tg.csv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE118523&format=file&file=GSE118523%5F20161109%5Fold%5Fwt%5Ftg%2Ecsv%2Egz'")
    os.system("wget -O GSE118523_20161109_young_wt_tg.csv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE118523&format=file&file=GSE118523%5F20161109%5Fyoung%5Fwt%5Ftg%2Ecsv%2Egz'")
    os.system("gunzip GSE118523_20161109_old_wt_tg.csv.gz")
    os.system("gunzip GSE118523_20161109_young_wt_tg.csv.gz")

def load_data(file_path):
    """Load dataset."""
    return pd.read_csv(file_path, index_col=0)

def perform_volcano_plot(data, output_prefix):
    """Create a volcano plot and save it."""
    # Ensure relevant columns are numeric
    data['log2fc'] = pd.to_numeric(data['log2fc'], errors='coerce')
    data['pval'] = pd.to_numeric(data['pval'], errors='coerce')

    # Define thresholds for significance
    log2fc_threshold = 1
    pval_threshold = 0.05

    # Create a column to classify genes based on significance
    data['significant'] = (data['pval'] < pval_threshold) & (data['log2fc'].abs() > log2fc_threshold)

    # Create the volcano plot
    plt.figure(figsize=(10, 6))
    plt.scatter(data['log2fc'], data['pval'], c='gray', alpha=0.6, label='Non-significant')
    plt.scatter(data.loc[data['significant'], 'log2fc'], 
                data.loc[data['significant'], 'pval'], 
                c='red', alpha=0.8, label='Significant')

    # Add lines for thresholds
    plt.axvline(x=log2fc_threshold, color='blue', linestyle='--', linewidth=1, label=f'Log2FC = {log2fc_threshold}')
    plt.axvline(x=-log2fc_threshold, color='blue', linestyle='--', linewidth=1)
    plt.axhline(y=pval_threshold, color='green', linestyle='--', linewidth=1, label=f'p-value = {pval_threshold}')

    # Customize plot
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('P-value')
    plt.title(f'Volcano Plot ({output_prefix})')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.savefig(f'{output_prefix}_volcano_plot.png')
    plt.close()

def perform_go_analysis(data, output_prefix):
    """Perform GO and pathway enrichment analysis and save results."""
    significant_genes = data.loc[data['significant'], 'gene_name'].dropna().tolist()
    if not significant_genes:
        print(f"No significant genes found for {output_prefix}.")
        return

    # Initialize g:Profiler object
    gp = GProfiler(return_dataframe=True)

    # Perform GO and pathway enrichment analysis
    go_results = gp.profile(
        organism='mmusculus',  # Mouse (Mus musculus)
        query=significant_genes,
        sources=['GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC']
    )

    # Filter results for significant enrichments
    go_results_filtered = go_results[go_results['p_value'] < 0.05]
    go_results_filtered.to_csv(f'{output_prefix}_go_pathway_results.csv', index=False)

    # Select top 10 GO terms by significance
    top_terms = go_results_filtered.sort_values('p_value').head(10)

    # Plot the results
    plt.figure(figsize=(10, 6))
    plt.barh(top_terms['name'], -np.log10(top_terms['p_value']), color='skyblue')
    plt.xlabel('-Log10(p-value)')
    plt.ylabel('GO Term')
    plt.title(f'Top 10 Enriched GO Terms ({output_prefix})')
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_go_terms.png')
    plt.close()

def main():
    # Download and prepare data
    download_and_prepare_data()

    # Analyze old mice dataset
    old_mice = load_data("GSE118523_20161109_old_wt_tg.csv")
    perform_volcano_plot(old_mice, "old_mice")
    perform_go_analysis(old_mice, "old_mice")

    # Analyze young mice dataset
    young_mice = load_data("GSE118523_20161109_young_wt_tg.csv")
    perform_volcano_plot(young_mice, "young_mice")
    perform_go_analysis(young_mice, "young_mice")

if __name__ == "__main__":
    main()

