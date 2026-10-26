import requests
import zipfile
import logging
import pandas as pd
from pathlib import Path

# ==========================================
# 1. LOGGING CONFIGURATION
# ==========================================
log_file_path = Path("data/logs/1_clean_metadata.log")
log_file_path.parent.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file_path, mode='w'),
        logging.StreamHandler()
    ]
)

def download_and_extract(url: str, base_dir: Path, zip_path: Path, expected_files: list) -> None:
    """Downloads and extracts the OSDR metadata zip file if it hasn't been downloaded yet."""
    # Check if files already exist
    if all(f.exists() for f in expected_files):
        logging.info("Metadata files already exist. Skipping download.")
        return

    if not zip_path.exists():
        logging.info("ZIP file not found. Downloading from OSDR...")
        try:
            response = requests.get(url, stream=True)
            response.raise_for_status()
            with open(zip_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            logging.info("Download complete.")
        except requests.exceptions.RequestException as e:
            logging.error(f"Failed to download the data: {e}")
            raise

    logging.info("Extracting ZIP file...")
    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(base_dir)
        logging.info("Extraction complete.")
    except zipfile.BadZipFile as e:
        logging.error(f"Failed to extract zip file: {e}")
        raise

def clean_metadata(input_file: Path, output_file: Path) -> None:
    """Loads, cleans, filters for Illumina Prep 3, and saves the metadata."""
    logging.info(f"Loading raw metadata from {input_file}...")
    try:
        df = pd.read_csv(input_file, sep='\t')
    except FileNotFoundError:
        logging.error(f"Input file not found: {input_file}. Did the ZIP extract correctly?")
        raise

    logging.info(f"Initial raw data shape: {df.shape}")

    cols_to_keep = [
        'Sample Name',
        'Characteristics[Biospecimen Type]',
        'Factor Value[Sample Location]',
        'Factor Value[Collection Date]',
        'Parameter Value[Prep Number]'
    ]

    # Subset and rename columns
    clean_df = df[cols_to_keep].copy()
    clean_df.columns = [
        'Sample_ID', 
        'Sample_Type', 
        'Location', 
        'Collection_Date', 
        'Prep_Number'
    ]

    logging.info("Sample Types distribution before filtering:")
    for key, value in clean_df["Sample_Type"].value_counts().items():
        logging.info(f"  {key}: {value}")

    # Standardize specific Sample ID
    clean_df['Sample_ID'] = clean_df['Sample_ID'].replace('Zymo_MMC_Log_11ng', 'Zymo_LD_MMC')

    # Filter for Prep 3 (Illumina short-read)
    logging.info("Filtering for Prep_Number == '3' (Illumina)...")
    clean_df = clean_df[
        clean_df['Prep_Number']
        .astype(str)
        .str.split(',')
        .apply(lambda prep_list: '3' in [item.strip() for item in prep_list])
    ].reset_index(drop=True)

    # Classify Sample_Category: distinguish true environmental cleanroom swabs from controls
    clean_df["Sample_Category"] = "Environmental"
    mask1 = clean_df["Sample_Type"].eq("Control")
    mask2 = (
        clean_df["Sample_Type"].eq("Environmental") &
        clean_df["Location"].str.contains("Control", case=False, na=False)
    )
    clean_df.loc[mask1 | mask2, "Sample_Category"] = "Control"

    logging.info("Sample Category breakdown:")
    for key, value in clean_df["Sample_Category"].value_counts().items():
        logging.info(f"  {key}: {value}")

    logging.info(f"Final cleaned data shape: {clean_df.shape}")

    # Ensure output directory exists and save
    output_file.parent.mkdir(parents=True, exist_ok=True)
    clean_df.to_csv(output_file, sep='\t', index=False)
    logging.info(f"Cleaned metadata successfully saved to {output_file}")

def main():
    # Define paths using pathlib for robust cross-operating-system compatibility
    base_dir = Path("data/archive")
    base_dir.mkdir(parents=True, exist_ok=True)
    
    zip_path = base_dir / "OSD-944_metadata_OSD-944-ISA.zip"
    file1 = base_dir / "a_OSD-944_metagenomic-sequencing_short-read-low-biomass-dna-sequencing_Illumina.txt"
    file2 = base_dir / "s_OSD-944.txt"
    output_file = Path("data/clean_metadata.tsv")
    
    download_url = "http://nasa-osdr.s3.amazonaws.com/OSD-944/version-1/metadata/OSD-944_metadata_OSD-944-ISA.zip"

    # Execute the pipeline
    download_and_extract(download_url, base_dir, zip_path, expected_files=[file1, file2])
    clean_metadata(file2, output_file)
    
    logging.info(f"Pipeline complete. Log saved to {log_file_path.absolute()}")

if __name__ == "__main__":
    main()

