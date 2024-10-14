import os
import subprocess
import argparse
import logging
import concurrent.futures

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def add_parser(subparsers):
    parser = subparsers.add_parser('sam2bam', help='Scan folder and convert sam to bam file.')
    parser.add_argument('-s', '--samtools', type=str, help='Path to samtools')
    parser.add_argument('-t', '--threads', type=int, default=8, help='Number of threads to use')

def run_command(command):
    tfile = os.path.basename(command.split(" ")[-1])
    logging.info(f"Starting sort: {tfile}")
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    logging.info(f"Finished sort: {tfile}")
    return result.stdout if result.returncode == 0 else result.stderr

def run(samtools, threads=8):
    commands = []
    for file in os.listdir("./"):
        if not file.lower().endswith(".sam"):
            continue
        sam_file = file
        bam_file = file.replace(".sam", ".sort.bam")
        commands.append(f"{samtools} view -h {sam_file} | {samtools} sort -o {bam_file}")
    
    if len(commands) == 0:
        logging.error("No SAM files found.")
        return
    
    logging.info(f"Converting {len(commands)} SAM files to BAM files.")
    logging.info(f"Using {threads} threads.")

    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(run_command, cmd): cmd for cmd in commands}
        for future in concurrent.futures.as_completed(futures):
            cmd = futures[future]
            try:
                result = future.result()
                if result:
                    logging.info(result)
            except Exception as exc:
                logging.error(f'{cmd} generated an exception: {exc}')
