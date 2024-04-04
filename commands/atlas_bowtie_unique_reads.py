import argparse
import subprocess
import os

def show_help():
    with open(os.path.expanduser("~/.atlaslogo"), 'r') as f:
        print(f.read())
    print("""
Usage: auto_unique_bam.py -i INPUT_SAM -o OUTPUT_BAM [-e ENV_NAME] [-h]

   -i   Input SAM file.
   -o   Output BAM file.
   -e   Anaconda environment name (default: base).
   -h   Show this help message.
""")

def main():
    parser = argparse.ArgumentParser(add_help=False)  # 禁用内置的帮助
    parser.add_argument("-i", "--input", type=str, required=True, help="Input SAM file.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output BAM file.")
    parser.add_argument("-e", "--env", type=str, default="base", help="Anaconda environment name (default: base).")
    parser.add_argument("-h", "--help", action="store_true", help="Show this help message.")
    
    args = parser.parse_args()
    
    if args.help:
        show_help()
        return

    # samtools命令
    cmd1 = ["samtools", "view", "-h", args.input]
    cmd2 = ["samtools", "sort", "-n", "-o", args.output + ".query.temp.sam"]
    subprocess.run(cmd1, stdout=subprocess.PIPE).stdout
    subprocess.run(cmd2)

    print("finish sorting")

    # unique_part命令（假定其为可执行命令）
    cmd_unique = ["unique_part", args.output + ".query.temp.sam", args.output + "_unique.sam"]
    subprocess.run(cmd_unique)

    print("finish unique")

    # 再次使用samtools
    cmd_sort = ["samtools", "sort", args.output + "_unique.sam", "-o", args.output + "_unique.sort.bam"]
    subprocess.run(cmd_sort)

    print("finish sorting")

    # 删除临时文件
    os.remove(args.output + ".query.temp.sam")
    os.remove(args.output + "_unique.sam")

    print("finish cleaning")

if __name__ == "__main__":
    main()
