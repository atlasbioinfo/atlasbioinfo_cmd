import argparse
import logging
import os
from pathlib import Path
import subprocess

def create_r_script(file_paths, sample_names, conditions, output_file):
    r_code = f'''
library(DESeq2)
library(tximport)
library(dplyr)

# Create file paths vector
files <- c({", ".join(f'"{f}"' for f in file_paths)})
names(files) <- c({", ".join(f'"{n}"' for n in sample_names)})

# Import RSEM files
txicase <- tximport(files, type = "rsem", txIn = TRUE, txOut = TRUE)

# Create sample metadata
samples <- data.frame(
    row.names = c({", ".join(f'"{n}"' for n in sample_names)}),
    condition = factor(c({", ".join(f'"{c}"' for c in conditions)}))
)

# Create DESeq2 dataset
ddsTxiRNA <- DESeqDataSetFromTximport(txicase, colData = samples, design = ~ condition)

# Calculate normalized FPKM values
normalizedRNA <- fpkm(ddsTxiRNA)

# Write results to CSV
write.csv(normalizedRNA, file = "{output_file}")
'''
    
    # 写入临时R脚本文件
    script_path = Path('temp_deseq2_script.R')
    with open(script_path, 'w') as f:
        f.write(r_code)
    return script_path

def run(input_dir, output_file, name_filter=None, rscript=None):
    # Get all RSEM files with optional name filter
    rsem_files = [f for f in os.listdir(input_dir) 
                  if f.endswith('.rsem') and 
                  (name_filter is None or name_filter in f)]
    if not rsem_files:
        logging.error("No RSEM files found in the input directory")
        return False
    
    # List all files with numbers
    logging.info("Found RSEM files:")
    for idx, file in enumerate(rsem_files, 1):
        logging.info(f"{idx}. {file}")
    
    # Get number of conditions
    while True:
        try:
            num_conditions = int(input("\nEnter the number of conditions: "))
            if num_conditions > 0:
                break
            logging.error("Please enter a positive number")
        except ValueError:
            logging.error("Please enter a valid number")
    
    # Get file assignments for each condition
    file_paths = []
    sample_names = []
    conditions = []
    
    for i in range(num_conditions):
        condition_name = f"condition{i+1}"
        file_indices = input(f"\nEnter file numbers for {condition_name} (comma-separated): ")
        
        # Process indices
        try:
            indices = [int(idx.strip()) for idx in file_indices.split(',')]
            for idx in indices:
                if 1 <= idx <= len(rsem_files):
                    file_name = rsem_files[idx-1]
                    file_paths.append(str(Path(input_dir) / file_name))
                    sample_names.append(file_name.split('.')[0])
                    conditions.append(condition_name)
                    logging.info(f"Added file '{file_name}' to condition '{condition_name}'")
                else:
                    logging.error(f"Invalid file number: {idx}")
                    return False
        except ValueError:
            logging.error("Invalid input format")
            return False
    
    try:
        logging.info("Creating R script for DESeq2 analysis...")
        script_path = create_r_script(file_paths, sample_names, conditions, output_file)
        
        logging.info("Starting DESeq2 analysis...")
        result = subprocess.run([rscript, str(script_path)], 
                              capture_output=True, 
                              text=True)
        
        if result.returncode != 0:
            logging.error(f"R script execution failed: {result.stderr}")
            return False
            
        # 删除临时R脚本
        script_path.unlink()
        
        logging.info(f"DESeq2 analysis completed successfully")
        logging.info(f"Results written to {output_file}")
        return True
        
    except Exception as e:
        logging.error(f"Error running DESeq2 analysis: {e}")
        return False

def add_parser(subparsers):
    parser = subparsers.add_parser('deseq2fpkm', help='Run DESeq2 analysis and calculate FPKM values')
    parser.add_argument('-i', '--input', type=str, default=".",
                      help='Input directory containing RSEM results')
    parser.add_argument('-o', '--output', type=str, required=True,
                      help='Output file path for normalized FPKM values')
    parser.add_argument('-f', '--filter', type=str,
                      help='Filter RSEM files by name (case-sensitive)')
    parser.set_defaults(func=run)
    return parser

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, default=".",
                      help='Input directory containing RSEM results')
    parser.add_argument('-o', '--output', type=str, required=True,
                      help='Output file path for normalized FPKM values')
    parser.add_argument('-f', '--filter', type=str,
                      help='Filter RSEM files by name (case-sensitive)')
    args = parser.parse_args()
    
    run(args.input, args.output, args.filter)

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
    main()
