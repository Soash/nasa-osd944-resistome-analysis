import pandas as pd
import requests
import re
import logging
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ==========================================
# 1. LOGGING CONFIGURATION
# ==========================================
log_file_path = Path("data/logs/5_resistome.log")
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
# 2. DATABASE FETCHING FUNCTIONS
# ==========================================
def fetch_kegg_amr_database(db_dir: Path) -> set:
    """Downloads and parses the official Antimicrobial Resistance KO IDs from KEGG."""
    logging.info("Fetching official ARG database from KEGG (ko01504)...")
    kegg_url = "http://rest.kegg.jp/link/ko/br:ko01504"
    kegg_file = db_dir / "kegg_amr_ko01504.txt"

    if not kegg_file.exists():
        try:
            response = requests.get(kegg_url, timeout=15)
            response.raise_for_status()
            with open(kegg_file, 'w') as f:
                f.write(response.text)
        except requests.exceptions.RequestException as e:
            logging.error(f"Failed to download KEGG database: {e}")
            return set()

    amr_ko_set = set()
    try:
        with open(kegg_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    parts = line.split('\t')
                    if len(parts) > 1:
                        ko_id = parts[1].replace('ko:', '')
                        amr_ko_set.add(ko_id)
        logging.info(f"Successfully loaded {len(amr_ko_set)} official Antibiotic Resistance KOs.")
    except Exception as e:
        logging.error(f"Error reading local KEGG database: {e}")

    return amr_ko_set

def fetch_bacmet_database(db_dir: Path) -> set:
    """Downloads and parses the BacMet biocide/metal resistance database."""
    logging.info("Fetching official BRG database from BacMet...")
    bacmet_url = "http://bacmet.biomedicine.gu.se/download/BacMet2_EXP.753.mapping.txt"
    bacmet_file = db_dir / "BacMet2_EXP.753.mapping.txt"

    if not bacmet_file.exists():
        try:
            response = requests.get(bacmet_url, timeout=15)
            response.raise_for_status()
            with open(bacmet_file, 'w') as f:
                f.write(response.text)
        except requests.exceptions.RequestException as e:
            logging.error(f"Failed to download BacMet database: {e}")
            return set()

    bacmet_genes = set()
    try:
        bacmet_df = pd.read_csv(bacmet_file, sep='\t')
        genes = bacmet_df['Gene_name'].dropna().astype(str).str.lower().tolist()

        all_genes = []
        for g in genes:
            all_genes.append(g)
            if '/' in g:
                all_genes.extend([x.strip() for x in g.split('/')])

        bacmet_genes = set(all_genes)
        logging.info(f"Successfully loaded {len(bacmet_genes)} official Biocide/Metal Resistance Genes.")
    except Exception as e:
        logging.error(f"Error parsing BacMet: {e}")

    return bacmet_genes

# ==========================================
# 3. HELPER EXPORT FUNCTIONS
# ==========================================
def matches_bacmet(function_text: str, bacmet_set: set) -> bool:
    """Checks if any distinct word in the string exists in the BacMet gene set."""
    words = set(re.findall(r'\b[a-z0-9]+\b', function_text))
    return not words.isdisjoint(bacmet_set)

def export_grouped_tsv(df: pd.DataFrame, output_tsv: Path, log_name: str) -> None:
    """Groups target genes by genus and exports count + KO ID/Function to a clean TSV."""
    if df.empty:
        logging.warning(f"No {log_name} genes found to export.")
        return

    temp_df = df.copy()
    temp_df['KO_combined'] = temp_df['KO_ID'].astype(str) + " (" + temp_df['KO_function_lower'].astype(str) + ")"

    grouped_df = temp_df.groupby('genus')['KO_combined'].agg(
        Unique_Count='nunique',
        Gene_List=lambda x: ', '.join(x.dropna().unique())
    ).reset_index()

    grouped_df = grouped_df.sort_values(by='Unique_Count', ascending=False)

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    grouped_df.to_csv(output_tsv, sep='\t', index=False, header=False)
    logging.info(f"Successfully exported {log_name} TSV to: {output_tsv}")

def export_summary_tsvs(location_burden: pd.Series, genus_burden: pd.Series, aci_mechanisms: pd.Series, result_dir: Path) -> None:
    """Exports the FULL statistical summary tables into clean TSV files."""
    loc_df = location_burden.reset_index()
    loc_df.columns = ['Location', 'Total_AMR_CPM']
    loc_df.to_csv(result_dir / "summary_location_burden.tsv", sep='\t', index=False)

    gen_df = genus_burden.reset_index()
    gen_df.columns = ['Genus', 'Total_AMR_CPM']
    gen_df.to_csv(result_dir / "summary_genera_burden.tsv", sep='\t', index=False)

    aci_df = aci_mechanisms.reset_index()
    aci_df.columns = ['Resistance_Mechanism', 'Gene_Count']
    aci_df.to_csv(result_dir / "summary_acinetobacter_mechanisms.tsv", sep='\t', index=False)
    
    logging.info("Successfully exported the full summary tables as TSV files.")

# ==========================================
# 4. ANALYSIS & PLOTTING FUNCTION
# ==========================================
def analyze_and_plot(env_amr_df: pd.DataFrame, fig_dir: Path, result_dir: Path) -> None:
    """Calculates summaries, logs global findings, exports TSVs, and generates plots."""
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Calculate Aggregates for logging
    total_amr_hits = len(env_amr_df)
    unique_genera = env_amr_df['genus'].nunique()
    unique_mechanisms = env_amr_df['KO_function'].nunique()

    # Calculate Data Summaries
    location_burden = env_amr_df.groupby('Location')['CPM'].sum().sort_values(ascending=False)
    genus_burden = env_amr_df.groupby('genus')['CPM'].sum().sort_values(ascending=False)
    
    acinetobacter_amr = env_amr_df[env_amr_df['genus'] == 'Acinetobacter']
    aci_mechanisms_full = acinetobacter_amr['KO_function'].value_counts()

    # Export 100% of the data to TSVs
    export_summary_tsvs(location_burden, genus_burden, aci_mechanisms_full, result_dir)

    # ---------------------------------------------------------
    # CLEAN LOGGING BLOCK
    # ---------------------------------------------------------
    logging.info("="*60)
    logging.info("           AMR ANALYSIS SUMMARY           ")
    logging.info("="*60)
    logging.info(f"    Total Environmental AMR Hits   : {total_amr_hits}")
    logging.info(f"    AMR-Carrying Bacterial Genera  : {unique_genera} unique genera")
    logging.info(f"    Unique AMR Mechanisms          : {unique_mechanisms} distinct functions")
    logging.info("    * Detailed breakdowns (Location, Genera, Acinetobacter) saved to TSV.")
    logging.info("="*60 + "\n")

    # Generate Plots
    logging.info(f"Generating and saving analytical plots to {fig_dir}...")
    sns.set_theme(style="whitegrid")
    plt.rcParams.update({'font.size': 12})
    
    genus_burden_top10 = genus_burden.head(10)
    aci_mechanisms_top5 = aci_mechanisms_full.head(5)

    # --- Figure 1: AMR Burden by Location ---
    plt.figure(figsize=(10, 6))
    sns.barplot(x=location_burden.values, y=location_burden.index, hue=location_burden.index, palette="Reds_r", legend=False)
    plt.title("AMR Burden by Cleanroom Location", fontsize=16, fontweight='bold')
    plt.xlabel("Normalized Abundance (Coverage Per Million)", fontsize=12)
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig(fig_dir / "amr_by_location.png", dpi=300)
    plt.close()

    # --- Figure 2: Top 10 AMR Carriers ---
    plt.figure(figsize=(10, 6))
    sns.barplot(x=genus_burden_top10.values, y=genus_burden_top10.index, hue=genus_burden_top10.index, palette="Blues_r", legend=False)
    plt.title("Top 10 AMR Carrying Bacteria", fontsize=16, fontweight='bold')
    plt.xlabel("Normalized Abundance (Coverage Per Million)", fontsize=12)
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig(fig_dir / "amr_by_genus.png", dpi=300)
    plt.close()

    # --- Figure 3: Acinetobacter Mechanisms ---
    plt.figure(figsize=(15, 5))
    sns.barplot(x=aci_mechanisms_top5.values, y=aci_mechanisms_top5.index, hue=aci_mechanisms_top5.index, palette="Purples_r", legend=False)
    plt.title("Top 5 Resistance Mechanisms in Acinetobacter", fontsize=16, fontweight='bold')
    plt.xlabel("Gene Hit Count", fontsize=12)
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig(fig_dir / "acinetobacter_mechanisms.png", dpi=300)
    plt.close()

# ==========================================
# 5. MAIN PIPELINE
# ==========================================
def main():
    db_dir = Path("data/database")
    db_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path("data/figures")
    result_dir = Path("data")
    
    matrix_file = result_dir / "master_functional_matrix.tsv"
    env_amr_file = result_dir / "env_amr_df.tsv"
    
    dual_tsv_file = result_dir / "total_dual_resistance_genes.tsv"
    amr_tsv_file = result_dir / "total_amr_genes.tsv"

    amr_ko_set = fetch_kegg_amr_database(db_dir)
    bacmet_genes = fetch_bacmet_database(db_dir)

    logging.info(f"Loading master matrix from {matrix_file}...")
    try:
        master_df = pd.read_csv(matrix_file, sep='\t')
    except FileNotFoundError:
        logging.error("Master matrix not found. Run previous pipeline steps first.")
        return

    logging.info("Cross-referencing genes against KEGG and BacMet databases...")
    master_df['KO_function_lower'] = master_df['KO_function'].fillna('').str.lower()
    
    master_df['is_BRG'] = master_df['KO_function_lower'].apply(lambda x: matches_bacmet(x, bacmet_genes))
    master_df["is_AMR"] = master_df["KO_ID"].isin(amr_ko_set)

    amr_resistance_df = master_df[master_df['is_AMR']]
    
    controls_to_remove = ['Zymo_LD_MMC', 'PC', 'Lambda_DCS']
    env_amr_df = amr_resistance_df[~amr_resistance_df['Sample_ID'].isin(controls_to_remove)].copy()
    
    logging.info(f"Filtered out mock communities. Final Environmental AMR gene hits: {len(env_amr_df)}\n")
    env_amr_df.to_csv(env_amr_file, sep='\t', index=False)

    dual_resistance_df = env_amr_df[env_amr_df['is_BRG'] & env_amr_df['is_AMR']]
    unique_dual_mechanisms = dual_resistance_df['KO_function'].nunique() if not dual_resistance_df.empty else 0

    logging.info("-"*60)
    logging.info("                 [ DUAL RESISTANCE ANALYSIS SUMMARY ]                ")
    logging.info("-" *60)
    logging.info(f"    Total Dual Resistance Gene Hits : {len(dual_resistance_df)}")
    
    if len(dual_resistance_df) > 0:
        unique_dual_counts = dual_resistance_df.groupby('genus')['KO_ID'].nunique()
        logging.info(f"    Dual-Resistant Genera Found     : {len(unique_dual_counts)} unique genera")
        logging.info(f"    Unique DUAL Resistance Mechs    : {unique_dual_mechanisms} distinct functions")
        logging.info("    * Detailed dual resistance lists saved to TSV.")
    else:
        logging.info("    0 genes identified carrying both AMR and BRG traits.")
    logging.info("-" *60 + "\n")

    export_grouped_tsv(env_amr_df, amr_tsv_file, "Total AMR")
    export_grouped_tsv(dual_resistance_df, dual_tsv_file, "Dual Resistance")

    analyze_and_plot(env_amr_df, fig_dir, result_dir)
    
    logging.info(f"Pipeline complete. All logs saved to {log_file_path.absolute()}")

if __name__ == "__main__":
    main()



