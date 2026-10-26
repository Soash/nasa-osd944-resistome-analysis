import requests
import logging
import pandas as pd
from pathlib import Path

# Set up the log file path
log_file_path = Path("data/logs/3_download_annotations.log")
log_file_path.parent.mkdir(parents=True, exist_ok=True)

# Configure logging to output timestamped status messages to BOTH terminal and a file
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file_path, mode='w'), # 'w' overwrites old logs; use 'a' to append
        logging.StreamHandler()
    ]
)

def download_resumable(url: str, dest_path: Path) -> bool:
    """Downloads a file from a URL, resuming if partially downloaded. Returns True if successful."""
    try:
        response = requests.head(url, timeout=10)
        response.raise_for_status()
        remote_size = int(response.headers.get('content-length', 0))
    except Exception as e:
        logging.error(f"Failed to fetch info for {url} | Error: {e}")
        return False

    local_size = 0
    if dest_path.exists():
        local_size = dest_path.stat().st_size

    if local_size == remote_size:
        logging.info(f"Already fully downloaded: {dest_path.name}")
        return True
    elif local_size > remote_size:
        logging.warning(f"Local file larger than remote. Restarting {dest_path.name}...")
        local_size = 0
        mode = 'wb'
    else:
        if local_size > 0:
            logging.info(f"Resuming {dest_path.name} from byte {local_size}...")
        else:
            logging.info(f"Starting {dest_path.name}...")
        mode = 'ab'

    headers = {'Range': f'bytes={local_size}-'} if local_size > 0 else {}

    try:
        with requests.get(url, headers=headers, stream=True, timeout=15) as r:
            r.raise_for_status()
            with open(dest_path, mode) as f:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
        logging.info(f"Successfully downloaded: {dest_path.name}")
        return True
    except Exception as e:
        logging.error(f"Error downloading {url} | Error: {e}")
        return False

def main():
    base_url = "http://nasa-osdr.s3-website-us-west-2.amazonaws.com/OSD-944/version-1/metagenomics-lowbiomass-shortread/"
    metadata_file = Path("data/clean_metadata.tsv")
    output_dir = Path("data/annotations")

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    logging.info(f"Reading sample IDs from {metadata_file}...")
    
    try:
        df = pd.read_csv(metadata_file, sep='\t')
    except FileNotFoundError:
        logging.error(f"Metadata file not found: {metadata_file}. Ensure previous pipeline steps ran successfully.")
        return

    # Safety check for column existence
    if 'Sample_ID' not in df.columns:
        logging.error("Column 'Sample_ID' not found in the cleaned metadata.")
        return

    target_samples = df['Sample_ID'].tolist()
    logging.info(f"Found {len(target_samples)} samples in metadata. Starting annotation downloads...")

    # Track successes and failures using lists
    successful_samples = []
    failed_samples = []

    # Loop through targets
    for sample in target_samples:
        filename = f"GLDS-775_metagenomics-lowbiomass-shortread_Prep3_{sample}-gene-coverage-annotation-and-tax.tsv"
        url = base_url + filename
        dest_path = output_dir / filename

        if download_resumable(url, dest_path):
            successful_samples.append(sample)
        else:
            failed_samples.append(sample)

    # --- FINAL SUMMARY BLOCK ---
    logging.info("=========================================")
    logging.info("         FINAL DOWNLOAD SUMMARY          ")
    logging.info("=========================================")
    logging.info(f"Total Expected Samples: {len(target_samples)}")
    logging.info(f"Successfully Acquired (New + Existing): {len(successful_samples)}")
    logging.info(f"Failed / Missing on Server: {len(failed_samples)}")
    logging.info("-----------------------------------------")
    
    if successful_samples:
        logging.info(f"SUCCESSFUL IDs: {', '.join(successful_samples)}")
    
    if failed_samples:
        logging.info(f"FAILED IDs: {', '.join(failed_samples)}")
        logging.warning("Note: Files like D1S7, D2S13, Buffer_AVE, and NTCs are known to be natively missing from the OSDR Prep3 server.")
    
    logging.info(f"Detailed log saved to: {log_file_path.absolute()}")

if __name__ == "__main__":
    main()


