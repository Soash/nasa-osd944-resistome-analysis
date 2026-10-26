import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import entropy
from scipy.spatial.distance import pdist, squareform
import seaborn as sns
import matplotlib.pyplot as plt
import logging
import warnings
from pathlib import Path

# Explicitly try to import scikit-bio (as it requires a separate installation)
try:
    from skbio.stats.ordination import pcoa
    from skbio.stats.distance import permanova
    from skbio import DistanceMatrix
    SKBIO_AVAILABLE = True
except ImportError:
    SKBIO_AVAILABLE = False

# ==========================================
# 1. LOGGING CONFIGURATION
# ==========================================
log_file_path = Path("data/logs/6_stat.log")
log_file_path.parent.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file_path, mode='w'),
        logging.StreamHandler()
    ]
)

# ==========================================
# 2. STATISTICAL FUNCTIONS
# ==========================================
def run_kruskal_wallis(env_amr_df: pd.DataFrame) -> pd.DataFrame:
    """Calculates total AMR load per sample and runs the Kruskal-Wallis H-test across locations."""
    logging.info("Calculating total AMR burden per sample...")
    sample_amr_load = env_amr_df.groupby(['Sample_ID', 'Location'])['CPM'].sum().reset_index()

    locations = sample_amr_load['Location'].unique()
    cpm_by_location = []

    for loc in locations:
        values = sample_amr_load[sample_amr_load['Location'] == loc]['CPM'].values
        if len(values) > 0:
            cpm_by_location.append(values)

    if len(cpm_by_location) > 1:
        stat, p_value = stats.kruskal(*cpm_by_location)
        
        logging.info("-"*50)
        logging.info("    [ KRUSKAL-WALLIS H-TEST RESULTS ]")
        logging.info("-" *50)
        logging.info(f"    H-statistic : {stat:.4f}")
        logging.info(f"    P-value     : {p_value:.4e}")
        
        if p_value < 0.05:
            logging.info("    Conclusion  : The difference in AMR burden between cleanroom locations is STATISTICALLY SIGNIFICANT.")
        else:
            logging.info("    Conclusion  : The difference is NOT statistically significant (p >= 0.05).")
        logging.info("-" *50)
    else:
        logging.warning("Not enough location groups to perform Kruskal-Wallis test.")

    return sample_amr_load

def analyze_alpha_diversity(env_amr_df: pd.DataFrame, sample_amr_load: pd.DataFrame, fig_dir: Path) -> pd.DataFrame:
    """Calculates Shannon Index (Alpha Diversity) and generates the corresponding plot."""
    logging.info("Calculating Alpha Diversity (Shannon Index)...")
    
    # Create Sample-by-Genus matrix
    genus_sample_matrix = env_amr_df.pivot_table(
        index='Sample_ID', 
        columns='genus', 
        values='CPM', 
        aggfunc='sum', 
        fill_value=0
    )

    # Function to calculate Shannon Index
    def shannon_index(row):
        proportions = row[row > 0] / row.sum()
        return entropy(proportions, base=np.e)

    # Apply to matrix
    alpha_diversity = genus_sample_matrix.apply(shannon_index, axis=1).reset_index()
    alpha_diversity.columns = ['Sample_ID', 'Shannon_Index']

    # Merge Location data
    alpha_diversity = alpha_diversity.merge(sample_amr_load[['Sample_ID', 'Location']], on='Sample_ID')

    # Generate Plot
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=alpha_diversity, x='Shannon_Index', y='Location', hue='Location', palette='viridis', fliersize=0, legend=False)
    sns.stripplot(data=alpha_diversity, x='Shannon_Index', y='Location', color='black', alpha=0.7, jitter=True)

    plt.title('Alpha Diversity (Shannon Index) of AMR Carriers by Location', fontsize=16, fontweight='bold')
    plt.xlabel('Shannon Diversity Index (Higher = More diverse resistance pool)', fontsize=12)
    plt.ylabel('')
    plt.tight_layout()
    
    output_img = fig_dir / "alpha_diversity.png"
    plt.savefig(output_img, dpi=300)
    plt.close()
    
    logging.info(f"Alpha Diversity plot saved to {output_img}")
    return genus_sample_matrix

def analyze_beta_diversity_and_permanova(genus_sample_matrix: pd.DataFrame, sample_amr_load: pd.DataFrame, fig_dir: Path):
    """Calculates Bray-Curtis distance, runs True PCoA, plots Beta Diversity, and executes PERMANOVA."""
    if not SKBIO_AVAILABLE:
        logging.error("The 'scikit-bio' library is not installed. Skipping Beta Diversity and PERMANOVA.")
        logging.info("To fix this, run: pip install scikit-bio")
        return

    logging.info("Calculating Bray-Curtis distance matrix for Beta Diversity...")
    bc_distance = pdist(genus_sample_matrix, metric='braycurtis')
    bc_distance_matrix = squareform(bc_distance)

    # Apply Square Root transformation to fix negative eigenvalues
    bc_distance_matrix_sqrt = np.sqrt(bc_distance_matrix)
    dm_sqrt = DistanceMatrix(bc_distance_matrix_sqrt, ids=genus_sample_matrix.index)
    dm_raw = DistanceMatrix(bc_distance_matrix, ids=genus_sample_matrix.index)

    # --- PCoA ---
    logging.info("Running Principal Coordinate Analysis (PCoA)...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        pcoa_results = pcoa(dm_sqrt, number_of_dimensions=2)

    pcoa_df = pcoa_results.samples[['PC1', 'PC2']].copy()
    pcoa_df['Sample_ID'] = pcoa_df.index
    pcoa_df = pcoa_df.merge(sample_amr_load[['Sample_ID', 'Location']], on='Sample_ID')

    pc1_var = pcoa_results.proportion_explained['PC1'] * 100
    pc2_var = pcoa_results.proportion_explained['PC2'] * 100

    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        data=pcoa_df, x='PC1', y='PC2', hue='Location', 
        s=150, palette='tab10', edgecolor='black', alpha=0.8
    )

    plt.title('Beta Diversity of Cleanroom AMR Communities', fontsize=16, fontweight='bold')
    plt.xlabel(f'Principal Coordinate 1 ({pc1_var:.1f}%)', fontsize=12)
    plt.ylabel(f'Principal Coordinate 2 ({pc2_var:.1f}%)', fontsize=12)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Location")
    plt.tight_layout()

    output_img = fig_dir / "beta_diversity_true_pcoa.png"
    plt.savefig(output_img, dpi=300)
    plt.close()
    logging.info(f"Beta Diversity (PCoA) plot saved to {output_img}")

    # --- PERMANOVA ---
    logging.info("Running PERMANOVA (999 permutations) to test Beta Diversity significance...")
    metadata = sample_amr_load.drop_duplicates(subset=['Sample_ID']).set_index('Sample_ID')
    
    np.random.seed(42) # Set seed for reproducibility
    permanova_results = permanova(dm_raw, metadata, column='Location', permutations=999)

    logging.info("-"*50)
    logging.info("    [ PERMANOVA RESULTS ]")
    logging.info("-" *50)
    logging.info(f"    Method      : {permanova_results['method name']}")
    logging.info(f"    Sample size : {permanova_results['sample size']}")
    logging.info(f"    Test stat   : {permanova_results['test statistic']:.4f}")
    logging.info(f"    P-value     : {permanova_results['p-value']:.4f}")
    
    if permanova_results['p-value'] < 0.05:
        logging.info("    Conclusion  : The microbial composition is SIGNIFICANTLY different across locations.")
    else:
        logging.info("    Conclusion  : No significant difference in composition across locations.")
    logging.info("-" *50 + "\n")


# ==========================================
# 3. MAIN PIPELINE
# ==========================================
def main():
    # Define directories and file paths
    fig_dir = Path("data/figures")
    fig_dir.mkdir(parents=True, exist_ok=True)
    
    env_amr_file = Path("data/env_amr_df.tsv")

    logging.info(f"Loading environmental AMR matrix from {env_amr_file}...")
    try:
        env_amr_df = pd.read_csv(env_amr_file, sep='\t')
    except FileNotFoundError:
        logging.error(f"Required file not found: {env_amr_file}. Run the previous resistome pipeline first.")
        return

    # Set the visual style for all plots
    sns.set_theme(style="whitegrid")
    plt.rcParams.update({'font.size': 12})

    # Step 1: Total AMR Load & Kruskal-Wallis Test
    sample_amr_load = run_kruskal_wallis(env_amr_df)

    # Step 2: Shannon Index (Alpha Diversity)
    genus_sample_matrix = analyze_alpha_diversity(env_amr_df, sample_amr_load, fig_dir)

    # Step 3: PCoA and PERMANOVA (Beta Diversity)
    analyze_beta_diversity_and_permanova(genus_sample_matrix, sample_amr_load, fig_dir)

    logging.info(f"Statistical pipeline complete. Logs saved to {log_file_path.absolute()}")

if __name__ == "__main__":
    main()