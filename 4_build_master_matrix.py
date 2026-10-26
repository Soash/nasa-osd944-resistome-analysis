import pandas as pd
import logging
from pathlib import Path

# ==========================================
# 1. LOGGING CONFIGURATION
# ==========================================
log_file_path = Path("data/logs/4_build_master_matrix.log")
log_file_path.parent.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file_path, mode='w'),
        logging.StreamHandler()
    ]
)

def build_master_matrix(annotation_dir: Path, metadata_file: Path, assaydata_file: Path, output_file: Path) -> None:
    """Reads annotation files, merges with metadata/assay data, calculates CPM, and saves the matrix."""
    
    # ==========================================
    # 2. LOAD METADATA & ASSAY DATA
    # ==========================================
    logging.info("Loading cleaned metadata and assay data...")
    try:
        metadata_df = pd.read_csv(metadata_file, sep='\t')
        assaydata_df = pd.read_csv(assaydata_file, sep='\t')
        logging.info(f"Successfully loaded metadata ({len(metadata_df)} samples) and assay data ({len(assaydata_df)} samples).")
    except FileNotFoundError as e:
        logging.error(f"Required input file missing: {e}")
        return

    # ==========================================
    # 3. PARSE ANNOTATION FILES
    # ==========================================
    annotation_files = list(annotation_dir.glob("*gene-coverage-annotation-and-tax.tsv"))
    
    if not annotation_files:
        logging.error(f"No annotation files found in {annotation_dir}. Check the download step.")
        return
        
    logging.info(f"Found {len(annotation_files)} annotation files. Processing and concatenating...")

    df_list = []
    target_cols = ['gene_ID', 'coverage', 'KO_ID', 'KO_function', 'genus', 'species']
    processed_samples = []

    for file in annotation_files:
        try:
            # Safely extract Sample_ID from the strict OSDR filename format
            sample_id = file.name.split('Prep3_')[1].split('-gene')[0]
            
            temp_df = pd.read_csv(file, sep='\t', usecols=target_cols)
            temp_df['Sample_ID'] = sample_id
            df_list.append(temp_df)
            processed_samples.append(sample_id)
            
        except IndexError:
            logging.error(f"Filename {file.name} does not match expected Prep3 format. Skipping.")
        except ValueError as e:
            logging.error(f"Missing expected columns in {file.name}: {e}. Skipping.")
        except Exception as e:
            logging.error(f"Unexpected error reading {file.name}: {e}")

    if not df_list:
        logging.error("No valid dataframes were created. Halting pipeline.")
        return

    # ==========================================
    # 4. CONCATENATE & FILTER
    # ==========================================
    logging.info("Concatenating all genes into a single master matrix...")
    master_df = pd.concat(df_list, ignore_index=True)
    logging.info(f"Raw master matrix shape: {master_df.shape}")

    logging.info("Filtering out genes without KO_ID annotations...")
    master_df = master_df.dropna(subset=['KO_ID'])
    logging.info(f"Filtered matrix shape (KO_ID present): {master_df.shape}")

    # ==========================================
    # 5. MERGE & NORMALIZE (CPM)
    # ==========================================
    logging.info("Merging functional data with physical cleanroom locations and assay depths...")
    
    # Inner join guarantees we only keep samples that exist in all 3 files
    meta_cols = ['Sample_ID', 'Location', 'Sample_Type']
    if 'Sample_Category' in metadata_df.columns:
        meta_cols.append('Sample_Category')
    final_df = master_df.merge(
        metadata_df[meta_cols], on='Sample_ID', how='inner'
    ).merge(
        assaydata_df[['Sample_ID', 'Total_Reads']], on='Sample_ID', how='inner'
    )

    logging.info("Calculating Coverage Per Million (CPM)...")
    final_df['CPM'] = (final_df['coverage'] / final_df['Total_Reads']) * 1000000

    logging.info(f"Final merged and normalized matrix shape: {final_df.shape}")
    
    # Verify if any samples were dropped during the inner join
    final_sample_count = final_df['Sample_ID'].nunique()
    if final_sample_count < len(processed_samples):
        logging.warning(f"Inner join dropped {len(processed_samples) - final_sample_count} samples. "
                        f"Check metadata/assaydata for missing Sample_IDs.")

    # ==========================================
    # 6. SAVE OUTPUT
    # ==========================================
    output_file.parent.mkdir(parents=True, exist_ok=True)
    final_df.to_csv(output_file, sep='\t', index=False)
    logging.info(f"Master functional matrix successfully saved to {output_file}")


def main():
    # Define paths using pathlib
    annotation_dir = Path("data/annotations")
    metadata_file = Path("data/clean_metadata.tsv")
    assaydata_file = Path("data/clean_assaydata.tsv")
    output_file = Path("data/master_functional_matrix.tsv")

    build_master_matrix(annotation_dir, metadata_file, assaydata_file, output_file)

if __name__ == "__main__":
    main()


    