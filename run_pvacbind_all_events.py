import subprocess
import os
import argparse
import logging
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json
import pandas as pd
from typing import Dict, List # <--- Added missing import for Dict and List

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - [%(funcName)s] - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# Global variable for HLA alleles map 
HLA_ALLELES_MAP = {}

def convert_allele(allele_str: str) -> str:
    """
    Convert a 5-character allele code (e.g., 'A0301') into 'HLA-A*03:01'.
    """
    if not isinstance(allele_str, str) or len(allele_str) != 5:
        logging.warning(f"Invalid allele string format: {allele_str}. Skipping conversion.")
        return allele_str # Return original if invalid
    locus = allele_str[0].upper()
    digits = allele_str[1:]
    if len(digits) != 4 or not digits.isdigit():
        logging.warning(f"Invalid allele digits in {allele_str}. Skipping conversion.")
        return allele_str
    group1 = digits[:2]
    group2 = digits[2:]
    return f"HLA-{locus}*{group1}:{group2}"

def convert_alleles_for_sample(allele_data: Dict[str, str]) -> str:
    """
    Convert a dictionary of sample HLA alleles into a comma-separated string required by pVACbind.
    E.g., {"A1": "A0301", "B1": "B3501"} -> "HLA-A*03:01,HLA-B*35:01"
    """
    converted_alleles = []
    for field, allele_str in allele_data.items():
        converted_allele = convert_allele(allele_str)
        # Only append if conversion was successful and it's not the original invalid string
        if converted_allele != allele_str or (converted_allele == allele_str and len(allele_str) == 5 and allele_str[1:].isdigit()):
            converted_alleles.append(converted_allele)
    return ",".join(converted_alleles)

def run_pvacbind_for_sample(fasta_file: Path, results_base_dir: Path, hla_alleles_map: Dict[str, str], 
                            pvacbind_docker_image: str, iedb_install_dir: str, 
                            binding_prediction_algorithms: List[str], threads_per_sample: int,
                            force_overwrite: bool = False):
    """
    Runs pvacbind for a single sample's FASTA file.

    Args:
        fasta_file (Path): Path to the input FASTA file with neoantigen peptides.
        results_base_dir (Path): Base directory for pVACbind results. A subdirectory will be created for each sample.
        hla_alleles_map (Dict[str, str]): Dictionary mapping sample_id to comma-separated HLA alleles string.
        pvacbind_docker_image (str): Name of the pVACbind Docker image.
        iedb_install_dir (str): Path to IEDB installation directory inside the Docker container.
        binding_prediction_algorithms (List[str]): List of binding prediction algorithms to use.
        threads_per_sample (int): Number of threads for pVACbind to use per sample.
        force_overwrite (bool): If True, overwrite existing pVACbind output for this sample.
    """
    # Assuming FASTA files are named like 'SAMPLE_ID_splicing_peptides.fasta'
    sample_name = fasta_file.stem.split('_splicing_peptides')[0] 
    sample_output_dir = results_base_dir / sample_name
    
    # Typical pVACbind output file to check for completion
    pvacbind_completion_file = sample_output_dir / f"{sample_name}.all_epitopes.tsv" 

    if pvacbind_completion_file.exists() and not force_overwrite:
        logging.info(f"Skipping {fasta_file.name}: pVACbind output already found in {sample_output_dir}. Use --force-overwrite to rerun.")
        return

    sample_output_dir.mkdir(parents=True, exist_ok=True)

    hla_alleles_str = hla_alleles_map.get(sample_name)
    if not hla_alleles_str:
        logging.warning(f"HLA alleles not found for sample {sample_name}. Skipping pVACbind.")
        return

    # Define paths for Docker volume mounts:
    # We need to mount the directory containing the fasta_file and the results directory.
    # The fasta_file itself needs to be referenced from within the Docker container's mount point.
    fasta_parent_dir = fasta_file.parent
    docker_fasta_mount_point = Path("/input_fasta")
    docker_results_mount_point = Path("/results")
    
    # Construct the input file path inside the Docker container
    docker_fasta_file = docker_fasta_mount_point / fasta_file.name

    algorithms_str = " ".join(binding_prediction_algorithms)

    docker_command = [
        "docker", "run", "--rm",
        "-v", f"{fasta_parent_dir}:{docker_fasta_mount_point}",
        "-v", f"{sample_output_dir}:{docker_results_mount_point}",
        "--user", f"{os.getuid()}:{os.getgid()}",
        pvacbind_docker_image,
        "pvacbind", "run",
        str(docker_fasta_file),
        sample_name,
        hla_alleles_str,
        algorithms_str,
        str(docker_results_mount_point),
        "--iedb-install-directory", iedb_install_dir,
        "-t", str(threads_per_sample)
    ]
    
    logging.info(f"Running pvacbind for sample {sample_name}...")
    logging.debug(f"pVACbind command: {' '.join(map(str, docker_command))}")
    
    try:
        subprocess.run(docker_command, check=True, capture_output=True, text=True)
        logging.info(f"pvacbind completed successfully for sample {sample_name}.")
        logging.debug(f"pvacbind stdout for {sample_name}:\n{result.stdout}")
    except subprocess.CalledProcessError as e:
        logging.error(f"pvacbind failed for sample {sample_name}. Error:\n{e.stderr}\nStdout:\n{e.stdout}")
        # Clean up partial output directory if error occurred
        if sample_output_dir.exists() and not list(sample_output_dir.iterdir()):
            try:
                sample_output_dir.rmdir()
            except OSError:
                pass # Directory might not be empty if some files were written
        raise
    except FileNotFoundError:
        logging.error("Docker command not found. Is Docker installed and in your PATH?")
        raise


def main():
    parser = argparse.ArgumentParser(
        description="Run pVACbind for all peptide FASTA files in a specified directory using Docker.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input-fasta-dir", type=Path, required=True,
                        help="Directory containing input peptide FASTA files (e.g., from process_spladder_output.py).")
    parser.add_argument("-o", "--output-dir", type=Path, required=True,
                        help="Root directory to save pVACbind results. A subdirectory will be created for each sample.")
    parser.add_argument("-a", "--hla-json-file", type=Path, required=True,
                        help="Path to a JSON file containing sample-HLA allele mappings.")
    parser.add_argument("--docker-image", type=str, default="griffithlab/pvactools",
                        help="Name of the pVACbind Docker image to use.")
    parser.add_argument("--iedb-install-dir", type=str, default="/opt/iedb",
                        help="Path to the IEDB installation directory *inside* the Docker container.")
    parser.add_argument("--algorithms", nargs='+', default=["NetMHCpan", "NetMHC", "DeepImmuno"],
                        help="List of binding prediction algorithms to use (e.g., NetMHCpan NetMHC DeepImmuno).")
    parser.add_argument("--workers", type=int, default=os.cpu_count() // 2, # Use half cores by default for Docker overhead
                        help="Number of parallel processes (Docker containers) to use.")
    parser.add_argument("--threads-per-sample", type=int, default=12,
                        help="Number of threads for pVACbind to use per sample/Docker container.")
    parser.add_argument("--force-overwrite", action="store_true",
                        help="Force rerun pVACbind even if output files already exist for a sample.")

    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    global HLA_ALLELES_MAP
    try:
        with open(args.hla_json_file, 'r') as f:
            raw_hla_dict = json.load(f)
        # Convert all HLA alleles once
        for sample_name, allele_data in raw_hla_dict.items():
            HLA_ALLELES_MAP[sample_name] = convert_alleles_for_sample(allele_data)
        logging.info(f"Loaded HLA alleles for {len(HLA_ALLELES_MAP)} samples.")
    except FileNotFoundError:
        logging.error(f"HLA JSON file not found: {args.hla_json_file}. Exiting.")
        return
    except json.JSONDecodeError as e:
        logging.error(f"Error parsing HLA JSON file {args.hla_json_file}: {e}. Exiting.")
        return
    except Exception as e:
        logging.error(f"An unexpected error occurred while loading HLA alleles: {e}. Exiting.")
        return

    # Discover FASTA files: process all found FASTA files
    fasta_files_to_process = sorted(list(args.input_fasta_dir.glob('*_splicing_peptides.fasta')))
    
    if not fasta_files_to_process:
        logging.warning("No FASTA files found to process. Please ensure your input directory contains files matching *_splicing_peptides.fasta. Exiting.")
        return

    logging.info(f"Found {len(fasta_files_to_process)} FASTA files to process.")
    logging.info(f"Using {args.workers} workers for pVACbind.")

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = []
        for fasta_file in fasta_files_to_process:
            futures.append(
                executor.submit(
                    run_pvacbind_for_sample, 
                    fasta_file, 
                    args.output_dir, 
                    HLA_ALLELES_MAP, 
                    args.docker_image, 
                    args.iedb_install_dir, 
                    args.algorithms, 
                    args.threads_per_sample,
                    args.force_overwrite # Pass the force_overwrite flag
                )
            )
        
        for future in futures:
            try:
                future.result()
            except Exception as e:
                logging.error(f"A pVACbind worker process failed: {e}", exc_info=True)

    logging.info("pVACbind pipeline finished.")

if __name__ == "__main__":
    main()