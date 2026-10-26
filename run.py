import subprocess
import sys
import logging
from pathlib import Path

# ==========================================
# 1. LOGGING CONFIGURATION
# ==========================================
# Setup logging to terminal to monitor the pipeline's overall progress
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - [PIPELINE] - %(levelname)s - %(message)s'
)

def run_command(command: list, error_message: str) -> None:
    """Helper function to execute shell commands safely."""
    try:
        # check=True ensures the pipeline halts immediately if a script crashes
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"{error_message} | Error Code: {e.returncode}")
        sys.exit(1)

def main():
    # 1. Define Paths
    venv_dir = Path("venv")
    
    # Handle executable paths depending on the Operating System
    if sys.platform == "win32":
        python_exe = venv_dir / "Scripts" / "python.exe"
        pip_exe = venv_dir / "Scripts" / "pip.exe"
    else:
        python_exe = venv_dir / "bin" / "python"
        pip_exe = venv_dir / "bin" / "pip"

    # 2. Create Virtual Environment
    if not venv_dir.exists():
        logging.info(f"Creating virtual environment at '{venv_dir.absolute()}'...")
        run_command(
            [sys.executable, "-m", "venv", str(venv_dir)], 
            "Failed to create virtual environment."
        )
    else:
        logging.info("Virtual environment already exists. Skipping creation.")

    # 3. Install Dependencies
    packages = ["pandas", "requests", "matplotlib", "seaborn", "scipy", "scikit-bio"]
    logging.info("Upgrading pip and installing required packages...")
    
    # Always a good idea to upgrade pip first
    run_command([str(python_exe), "-m", "pip", "install", "--upgrade", "pip"], "Failed to upgrade pip.")
    
    # Install the actual requirements
    run_command([str(pip_exe), "install"] + packages, "Failed to install dependencies.")

    # 4. Execute the Analytical Pipeline
    scripts = [
        "1_clean_metadata.py",
        "2_clean_assaydata.py",
        "3_download_annotations.py",
        "4_build_master_matrix.py",
        "5_resistome.py",
        "6_stat.py"
    ]

    logging.info("========================================")
    logging.info("      STARTING ANALYTICAL PIPELINE      ")
    logging.info("========================================")

    for script in scripts:
        script_path = Path(script)
        if not script_path.exists():
            logging.error(f"Script missing: {script_path}. Halting pipeline.")
            sys.exit(1)

        logging.info(f"Executing: {script} ...")
        # Running the script using the virtual environment's Python
        run_command([str(python_exe), str(script_path)], f"Pipeline failed during {script}")
        logging.info(f"Successfully completed: {script}")

    logging.info("========================================")
    logging.info("      PIPELINE COMPLETED SUCCESSFULLY   ")
    logging.info("========================================")

if __name__ == "__main__":
    main()


    