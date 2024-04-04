#!/Users/hyu/miniforge3/bin/python

import argparse
import os,tqdm

def check_big_file(dir_path, size_threshold):
    target_file_count=0
    max_file_size=-1
    for root, dirs, files in tqdm.tqdm(os.walk(dir_path)):
        for file in files:
            file_path = os.path.join(root, file)
            file_size = os.path.getsize(file_path)
            if file_size > size_threshold*1024*1024:
                target_file_count+=1
                if file_size > max_file_size:
                    max_file_size = file_size
                # print(f"{file_path} is {file_size/1024/1024:.2f}MB", end="\r")
                if compress:
                    os.system(f"pigz {file_path}")
                    print(f"{file_path} is compressed")
    print(f"Total {target_file_count} files are larger than {size_threshold}MB")
    print("Done")
                


if __name__ == "__main__":

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
    This script is used to check the size of all files in a directory.
    If the file size is larger than 100MB, it will be printed out.
    -z to compress the output file.
    .'''.format(logo)

    parser = argparse.ArgumentParser(description=description_text, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-s', '--size', type=int, default=100, help='file size threshold in MB, default is 100MB')
    #compress or not, default is False
    parser.add_argument('-z', '--compress', action='store_true', help='compress the file')
    
    args = parser.parse_args()
    size_threshold = args.size
    compress = args.compress
    # print(f"size_threshold: {size_threshold}MB")
    # print(f"compress: {compress}")
    dir_path = os.getcwd()
    check_big_file(dir_path, size_threshold)
    