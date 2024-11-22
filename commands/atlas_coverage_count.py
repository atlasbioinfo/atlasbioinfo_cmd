#!/Users/hyu/Intel_miniforge3/envs/PYSAM/bin/python
# -*- coding: utf-8 -*-
import argparse
from pathlib import Path
import logging

def checkread(read): #reverse strand
    # if not read.is_paired:
    if (read.is_reverse):
        return False
    # if read.is_read1 and read.mate_is_reverse:
    #     return True
    # if read.is_read2 and read.mate_is_reverse:
    #     return True
    
    return True

def combine4Arr(arrs):
    cArr=[]
    for i in range(len(arrs[0])):
        cArr.append(str(arrs[0][i]+arrs[1][i]+arrs[2][i]+arrs[3][i]))
    return ",".join(cArr)

def main(input_file, output, reads_all, fasta):
    import os,re,sys,pysam,time
    import numpy as np
    
    if (fasta is not None):
        seq={}
        with open(fasta) as f:
            fcont=f.read().split(">")
            for cont in fcont[1:]:
                tmp=cont.split("\n")
                tname=re.split("\s+",tmp[0])[0]
                # tname=tmp[0].split(" ")[0]
                tseq="".join(tmp[1:])
                seq[tname]=tseq.upper()

    beg=time.time()
    bamfile=pysam.AlignmentFile(input_file, "rb")

    discardReverse=0
    calculate_reads=0
    # count=0
    # limit=100000
    with open(output,"w") as out:
        for i in range(bamfile.nreferences):
            tname=bamfile.get_reference_name(i)
            tlen=bamfile.get_reference_length(tname)
            if (bamfile.count(tname)==0):
                continue
            genePos=[0 for i in range(int(tlen))]
            geneCov=np.array([0 for i in range(int(tlen))])
            for read in bamfile.fetch(tname):
                if (read.is_unmapped):
                    continue
                calculate_reads+=1
                if (reads_all):
                    genePos[read.reference_start]+=1
                    geneCov[read.reference_start:read.reference_end+1]+=1
                    continue
                if (checkread(read)):
                    # count+=1
                    # if (count>limit):
                    #     break
                    discardReverse+=1
                    genePos[read.reference_start]+=1
                    geneCov[read.reference_start:read.reference_end+1]+=1

            genePos=[str(i) for i in genePos]
            geneCov=[str(i) for i in geneCov]

            if (fasta is not None):
                out.write("\t".join([
                    tname,
                    seq[tname],
                    ",".join(genePos),
                    ",".join(geneCov)
                ])+"\n")
            else:    
                out.write("\t".join([
                    tname,
                    ",".join(genePos),
                    ",".join(geneCov)
                ])+"\n")
    
    if (reads_all):
        print("Total reads: "+str(calculate_reads))
    else:
        print("Total reads: "+str(calculate_reads))
        print("Calculated (discard reverse reads): "+str(discardReverse))
    
    print("Time elapsed: "+str(time.time()-beg)+"s")

def run(pysam_python, args):
    import subprocess
    if pysam_python is None:
        print("Error: pysam_python path is None. Please make sure the pysam environment is properly set up.")
        return

    script_path = Path(__file__).resolve()
    command = [pysam_python, str(script_path)]  # Start with base command
    command.append(args.input_bam)

    if args.output:
        command.extend(['-o', args.output])
    else:
        default_output = args.input_bam + '.AtlasCovRT'
        command.extend(['-o', default_output])
    
    if args.fasta:
        command.extend(['-f', args.fasta])
    
    if args.all:
        command.extend(['-a', str(args.all)])
    # print(" ".join(command))
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running atlas coverage count: {e}")
    except FileNotFoundError:
        print(f"Error: Could not find the pysam Python interpreter at {pysam_python}")

def add_parser(subparsers):
    parser = subparsers.add_parser('cov_rt_count', help='Calculate coverage and position counts for each reference in a BAM file')
    parser.add_argument('input_bam', type=str, help='Path to the input indexed BAM file.')
    parser.add_argument('-f', '--fasta', type=str, default=None, help='Path to the reference fasta file, default is None')
    parser.add_argument('-o', '--output', type=str, default=None, help='output file name, default is input_file.AtlasCovRT')
    parser.add_argument('-a', '--all', type=bool, default=False, help='output all reads, include the antisense stranded reads, default is False')
    parser.set_defaults(func=run)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
    parser = argparse.ArgumentParser(description='Calculate coverage and position counts for each reference in a BAM file')
    parser.add_argument('input_bam', type=str, help='Path to the input indexed BAM file.')
    parser.add_argument('-f', '--fasta', type=str, default=None, help='Path to the reference fasta file, default is None')
    parser.add_argument('-o', '--output', type=str, default=None, help='output file name, default is input_file.AtlasCovRT')
    parser.add_argument('-a', '--all', type=bool, default=False, help='output all reads, include the antisense stranded reads, default is False')
    args = parser.parse_args()

    if args.output is None:
        args.output = args.input_bam + '.AtlasCovRT'
    
    main(args.input_bam, args.output, args.all, args.fasta)
