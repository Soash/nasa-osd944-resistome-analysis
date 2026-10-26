import logging
import pandas as pd
from pathlib import Path

# ==========================================
# 1. LOGGING CONFIGURATION
# ==========================================
log_file_path = Path("data/logs/2_clean_assaydata.log")
log_file_path.parent.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file_path, mode='w'),
        logging.StreamHandler()
    ]
)

def clean_assay_data(input_file: Path, output_file: Path) -> None:
    """Loads raw assay data, filters for specific QA and depth metrics, and saves."""
    logging.info(f"Loading raw assay data from {input_file}...")
    
    try:
        df = pd.read_csv(input_file, sep='\t')
    except FileNotFoundError:
        logging.error(f"Input file not found: {input_file}. Ensure the OSDR archive was downloaded and extracted.")
        raise

    logging.info(f"Initial raw assay data shape: {df.shape}")

    # Define columns to extract based on OSDR metadata structure
    cols_to_keep = [
        'Sample Name',
        'Parameter Value[QA Score]',
        'Parameter Value[Read Depth]',
        'Parameter Value[Human Reads Removed]',
        'Parameter Value[Contaminant Reads Removed]'
    ]

    # Safety Check: Verify all expected columns exist to prevent KeyError
    missing_cols = [col for col in cols_to_keep if col not in df.columns]
    if missing_cols:
        logging.error(f"Missing expected columns in input file: {missing_cols}")
        raise KeyError(f"Missing columns: {missing_cols}")

    # Subset the dataframe
    clean_df = df[cols_to_keep].copy()

    # Rename columns to standardized pipeline format
    clean_df.columns = [
        'Sample_ID', 
        'DNA_Concentration_ng_ul', 
        'Total_Reads', 
        'Human_Reads_Removed_Pct', 
        'Blank_Contaminant_Removed_Pct'
    ]

    # Standardize the Zymo Log-Distributed Mock Community name
    clean_df['Sample_ID'] = clean_df['Sample_ID'].replace('Zymo_MMC_Log_11ng', 'Zymo_LD_MMC')

    logging.info(f"Final cleaned assay data shape: {clean_df.shape}")

    # Ensure output directory exists and save
    output_file.parent.mkdir(parents=True, exist_ok=True)
    clean_df.to_csv(output_file, sep='\t', index=False)
    logging.info(f"Cleaned assay data successfully saved to {output_file}")

def main():
    # Define paths using pathlib for robust cross-OS compatibility
    input_file = Path("data/archive/a_OSD-944_metagenomic-sequencing_short-read-low-biomass-dna-sequencing_Illumina.txt")
    output_file = Path("data/clean_assaydata.tsv")
    
    clean_assay_data(input_file, output_file)
    
    logging.info(f"Pipeline complete. Log saved to {log_file_path.absolute()}")

if __name__ == "__main__":
    main()