import os
import subprocess
import argparse
import logging
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - [%(funcName)s] - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def run_spladder_build(bam_file: Path, output_dir: Path, gtf_file: Path) -> None:
    """
    Creates splice graphs for a single sample using SplAdder.

    Args:
        bam_file (Path): Path to the BAM file containing the RNA-seq reads.
        output_dir (Path): Path to the root directory where the splice graphs will be saved.
                           A subdirectory named after the sample will be created here.
        gtf_file (Path): Path to the GTF file containing the gene annotations.
        force_overwrite (bool): If True, overwrite existing SplAdder output for this sample.
    """
    # Extract sample name
    sample_name = bam_file.stem.split('.')[0] 
    
    sample_output_dir = output_dir / sample_name

    sample_output_dir.mkdir(parents=True, exist_ok=True)

    logging.info(f"Creating splice graphs for {bam_file.name} in {sample_output_dir}...")
    cmd = [
        "spladder", "build",
        "-b", str(bam_file),
        "-a", str(gtf_file),
        "-o", str(sample_output_dir),
        "-M", "single",
        "--no-extract-ase",
        "--sparse-bam"
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logging.info(f"Successfully created splice graphs for {bam_file.name}")
    except subprocess.CalledProcessError as e:
        logging.error(f"SplAdder build failed for {bam_file.name}. Error:\n{e.stderr}\nStdout:\n{e.stdout}")
        # Clean up partial output directory if error occurred, to ensure rerun next time
        if sample_output_dir.exists() and not list(sample_output_dir.iterdir()):
            try:
                sample_output_dir.rmdir()
            except OSError:
                pass # Directory might not be empty if some files were written
        raise
    except FileNotFoundError:
        logging.error("SplAdder command not found. Is SplAdder installed and in your PATH?")
        raise

def main():
    parser = argparse.ArgumentParser(
        description="Run SplAdder 'build' command for all BAM files in a specified directory.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-b", "--bam-dir", type=Path, required=True,
                        help="Directory containing BAM files to process.")
    parser.add_argument("-o", "--output-dir", type=Path, required=True,
                        help="Root directory to save SplAdder output. A subdirectory will be created for each sample.")
    parser.add_argument("-g", "--gtf-file", type=Path, required=True,
                        help="Path to the GTF file containing the gene annotations.")
    parser.add_argument("--workers", type=int, default=os.cpu_count(),
                        help="Number of parallel processes to use.")

    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Add to list all the BAM files, ensure they have a .bai index for robust processing
    bam_files_to_process = []
    for bam_file in sorted(list(args.bam_dir.glob('*.bam'))):
        if bam_file.name.endswith('.out.bam'): 
            logging.debug(f"Skipping {bam_file.name}: Looks like an unmapped output BAM (ends with .out.bam).")
            continue
        if not bam_file.with_suffix('.bam.bai').exists():
            logging.warning(f"Skipping {bam_file.name}: No corresponding .bai index file found. BAM files require an index.")
            continue
        bam_files_to_process.append(bam_file)
        
    logging.info(f"Found {len(bam_files_to_process)} BAM files to process.")

    if not bam_files_to_process:
        logging.warning("No valid BAM files found to process. Please ensure BAM files have .bai indices. Exiting.")
        return

    logging.info(f"Using {args.workers} workers.")

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = []
        for bam_file in bam_files_to_process:
            futures.append(
                executor.submit(run_spladder_build, bam_file, args.output_dir, args.gtf_file)
            )
        
        for future in futures:
            try:
                future.result()
            except Exception as e:
                logging.error(f"A SplAdder worker process failed: {e}", exc_info=True)

    logging.info("SplAdder 'build' pipeline finished.")

if __name__ == "__main__":
    main()