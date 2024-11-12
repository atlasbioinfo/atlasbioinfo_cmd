import os,pysam,argparse
import subprocess
from pathlib import Path
import logging

def is_bam_file(filepath):
    try:
        with pysam.AlignmentFile(filepath, "rb") as bam_file:
            return True
    except ValueError:
        return False

def is_bam_indexed(filepath):
    try:
        bam_file = pysam.AlignmentFile(filepath, "rb")
        bam_index = pysam.IndexedReads(bam_file)
        bam_index.build()
        return True
    except FileNotFoundError:
        
        return False

def check_bam_file(input_bam):
    if not os.path.exists(input_bam):
        print("File not found: ", input_bam)
        return False

    try:
        with pysam.AlignmentFile(input_bam, "rb") as bam_file:
            print(input_bam, " is a valid BAM file")
    except ValueError:
        print(input_bam, " is not a valid BAM file")
        return False
    
    if not os.path.exists(input_bam+".bai"):
        print(input_bam, " is not indexed, creating index...")
        try:
            pysam.index(input_bam)
            print("Index created successfully")
        except Exception as e:
            print(f"Error creating index: {e}")
            return False
    else:
        print(input_bam, " is indexed")
    return True

def check_read(read, count_strand, dna):
    if count_strand == "B":  # Count both strands
        return True
        
    if dna:
        if read.is_paired:
            is_forward = (read.is_read1 and not read.is_reverse) or (read.is_read2 and read.is_reverse)
        else:
            #可能有问题
            is_forward = not read.is_reverse
    else:  
        if read.is_paired:
            is_forward = (read.is_read1 and read.is_reverse) or (read.is_read2 and not read.is_reverse)
        else:
            #可能有问题
            is_forward = not read.is_reverse
    
    if count_strand == "F":  # Forward strand
        return is_forward
    elif count_strand == "R":  # Reverse strand
        return not is_forward
    
    return False

def run_test(input_bam, count_strand="F", output_file=None, dna=False):
    if not check_bam_file(input_bam):
        return False
    
    if output_file is None:
        output_file = input_bam + ".test.count"
    
    bamfile = pysam.AlignmentFile(input_bam, "rb")
    
    # 统计变量初始化
    all_reads = 0
    unmapped_count = 0
    counted_reads = 0
    paired_reads = 0
    single_reads = 0
    test_limit = 10000  # 测试读取的reads数量限制
    
    print("\nTest Mode - Sampling first", test_limit, "reads...")
    
    for read in bamfile.fetch(until_eof=True):
        if all_reads >= test_limit:
            break
            
        all_reads += 1
        
        # 统计配对和单端reads
        if read.is_paired:
            paired_reads += 1
        else:
            single_reads += 1
            
        # 统计未比对reads
        if read.is_unmapped:
            unmapped_count += 1
            continue
            
        # 根据链特异性统计
        if check_read(read, count_strand, dna):
            counted_reads += 1
    
    # 计算统计结果
    mapped_rate = round((all_reads - unmapped_count) / all_reads * 100, 2)
    counted_rate = round(counted_reads / all_reads * 100, 2)
    paired_rate = round(paired_reads / all_reads * 100, 2)
    
    # 输出统计结果
    print("\nTest Results:")
    print(f"Total reads sampled: {all_reads}")
    print(f"Mapped rate: {mapped_rate}%")
    print(f"Counted reads: {counted_reads} ({counted_rate}%)")
    print(f"Library type:")
    print(f"  - Paired-end reads: {paired_reads} ({paired_rate}%)")
    print(f"  - Single-end reads: {single_reads} ({100-paired_rate}%)")
    print(f"Strand counting mode: {count_strand}")
    
    bamfile.close()
    return True

def run(input_bam, count_strand="F", output_file=None, dna=False):
    logging.info("Running count mode...")
    if not check_bam_file(input_bam):
        return False
    

    if output_file is None:
        output_file = input_bam + ".count"
    
    bamfile=pysam.AlignmentFile(input_bam, "rb")
    
    all_reads=0
    unmapped_count=0
    counted_reads=0
    gene_count=0
    with open(output_file,"w") as out:
        out.write("Gene\tLength\tCount\n")
        for i in range(bamfile.nreferences):
            tname=bamfile.get_reference_name(i)
            tlen=str(bamfile.get_reference_length(tname))
            if (bamfile.count(tname)==0):
                continue
            gene_count+=1
            tcount=0
            # if (all_reads>limit):
            #     break
            for read in bamfile.fetch(tname):
                all_reads+=1
                if (read.is_unmapped):
                    unmapped_count+=1
                    continue
                if check_read(read,count_strand,dna):
                    counted_reads+=1
                    tcount+=1
            out.write("\t".join([tname,tlen,str(tcount)])+"\n")

    print("Total reads: ",all_reads)
    print("Gene count: ",gene_count)
    print("Mapped rate: ",str(round((all_reads-unmapped_count)/all_reads,4)*100)+"%")
    print("Counted reads: ",counted_reads," Rate: ",str(round(counted_reads/all_reads,4)*100)+"%")

def run_count_rnaseq_with_pysam(pysam_python, args):
    if pysam_python is None:
        print("Error: pysam_python path is None. Please make sure the pysam environment is properly set up.")
        return

    script_path = Path(__file__).resolve()
    command = [
        pysam_python,
        str(script_path),
        args.input_bam,
    ]

    if hasattr(args, 'dna') and args.dna:
        command.append('-d')
    
    if hasattr(args, 'count_strand'):
        command.extend(['-s', args.count_strand])
    
    if hasattr(args, 'output') and args.output:
        command.extend(['-o', args.output])
    
    if hasattr(args, 'test') and args.test:
        command.append('-t')

    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running count_rnaseq: {e}")
    except FileNotFoundError:
        print(f"Error: Could not find the pysam Python interpreter at {pysam_python}")

def add_parser(subparsers):
    parser = subparsers.add_parser('count_rnaseq', help='Count reads in RNA-seq BAM file')
    parser.add_argument('input_bam', type=str, help='Path to the input indexed BAM file.')
    parser.add_argument('-s', '--count_strand', type=str, default="F", choices=['F', 'R', 'B'],
                        help='Strand to count: F (forward), R (reverse), B (both). Default: F')
    parser.add_argument('-o', '--output', type=str, help='Path to the output file. If not specified, uses input_bam + ".count"')
    parser.add_argument('-t', '--test', action='store_true',
                       help='Test mode: randomly sample 10000 reads for assessment')
    parser.add_argument('-d', '--dna', action='store_true',
                       help='Indicate if reads are mapped to DNA sequence')

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
    parser = argparse.ArgumentParser(description='Count reads in RNA-seq BAM file')
    parser.add_argument('input_bam', type=str, help='Path to the input indexed BAM file.')
    parser.add_argument('-s', '--count_strand', type=str, default="F", choices=['F', 'R', 'B'],
                        help='Strand to count: F (forward), R (reverse), B (both). Default: F')
    parser.add_argument('-o', '--output', type=str, help='Path to the output file. If not specified, uses input_bam + ".count"')
    parser.add_argument('-t', '--test', action='store_true',
                       help='Test mode: randomly sample 10000 reads for assessment')
    parser.add_argument('-d', '--dna', action='store_true',
                       help='Indicate if reads are mapped to DNA sequence')
    args = parser.parse_args()

    if args.output is None:
        args.output = args.input_bam + ".count"

    if args.test:
        logging.info("Running test mode...")
        run_test(args.input_bam, args.count_strand, args.output, args.dna)
    else:
        run(args.input_bam, args.count_strand, args.output, args.dna)
