#!/Users/hyu/Intel_miniforge3/envs/PYSAM/bin/python
# -*- coding: utf-8 -*-
import os
import subprocess, argparse

def add_parser(subparsers):
    parser = subparsers.add_parser('bamindex', help='Index BAM files')
    parser.set_defaults(func=bamindex)

def bamindex(args):
    run(args.samtools_path)

def run(samtools_path):
    count = 0
    countHasNoIndex = 0
    for file in os.listdir("./"):
        if file.endswith(".bam"):
            if not os.path.exists(file + ".bai"):
                countHasNoIndex += 1
                subprocess.call("nohup " + samtools_path + " index " + file + " &", shell=True)
            else:
                count += 1

    print("Total bam files: " + str(count + countHasNoIndex))
    print("Bam files with index: " + str(count))
    print("Bam files without index: " + str(countHasNoIndex))
    print("Indexing...Please check the backend process.")
