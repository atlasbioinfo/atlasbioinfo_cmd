#!/Users/hyu/Intel_miniforge3/envs/PYSAM/bin/python
# -*- coding: utf-8 -*-

import os,argparse
import time

def countdown(seconds):
    while seconds > 0:
        time.sleep(1)
        seconds -= 1
        print("Start in ",seconds," seconds",end="\r")
        
def main():
    print("This script is a pre-step for DESEQ2 to convert all counts to rsem format files, please make sure your data is RNA-Seq or Polysome data.")
    countdown(3)
    genelist={}
    file_count=0
    for file in os.listdir("./"):
        if (not file.endswith(".count")):
            continue
        file_count+=1
        with open(file,"r")as f:
            f.readline()
            for line in f:
                tmp=line.strip().split("\t")
                if (tmp[0] not in genelist):
                    genelist[tmp[0]]=0
                genelist[tmp[0]]+=1

    target_gene=set()
    for gene in genelist:
        if (genelist[gene]==file_count):
            target_gene.add(gene)

    for file in os.listdir("./"):
        if (not file.endswith(".count")):
            continue
        tname=file.split(".")[0]
        with open(tname+".rsem","w")as out:
            with open(file,"r")as f:
                f.readline()
                out.write("transcript_id\tgene_id\tlength\teffective_length\texpected_count\tTPM\tFPKM\tIsoPct\n")
                for line in f:
                    tmp=line.strip().split("\t")
                    if (tmp[0] not in target_gene):
                        continue
                    out.write("\t".join([
                        tmp[0],tmp[0],
                        tmp[1],tmp[1],
                        tmp[2],"0","0","0"
                    ])+"\n")

def add_parser(subparsers):
    parser = subparsers.add_parser('count2rsem', help='Convert all counts to rsem format files')

            
if __name__=="__main__":

    logo='''      
          _   _             ____  _       _        __      
     /\  | | | |           |  _ \(_)     (_)      / _|     
    /  \ | |_| | __ _ ___  | |_) |_  ___  _ _ __ | |_ ___  
   / /\ \| __| |/ _` / __| |  _ <| |/ _ \| | '_ \|  _/ _ \ 
  / ____ \ |_| | (_| \__ \ | |_) | | (_) | | | | | || (_) |
 /_/    \_\__|_|\__,_|___/ |____/|_|\___/|_|_| |_|_| \___/  

        `-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"
        `=`,'=/     `=`,'=/     `=`,'=/     `=`,'=/
            y==/        y==/        y==/        y==/
        ,=,-<=`.    ,=,-<=`.    ,=,-<=`.    ,=,-<=`.
        ,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_
                
    '''


    description_text = '''{} 
        "atlas_count_to_rsem" is a tool to convert all counts to rsem format files, please make sure your data is RNA-Seq or Polysome data.'''.format(logo)

    parser = argparse.ArgumentParser(description=description_text, formatter_class=argparse.RawTextHelpFormatter)
    args = parser.parse_args()
    
    
    main()