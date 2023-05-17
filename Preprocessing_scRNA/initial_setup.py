import os
import sys
import argparse


def main(args):
    # first check if the folder already exists
    dir_path = os.path.join(args.work_dir, "id" + args.id)
    if os.path.exists(dir_path):
        sys.exist("The directory for id", args.id, " already exists.")
    else:
        os.mkdir(dir_path)
        readme = open(os.path.join(dir_path, "readme.md"), "w")
        print("#Documentation", file=readme)

def parse_args():
    parser = argparse.ArgumentParser(description='Initial setup for processing new scRNAseq data')
    parser.add_argument('--id', required=True, help='Input id as a number. For example: 1, 2, 3, 4, etc...')
    parser.add_argument('--work_dir', required=True, help='Input working directory where the folder is going to be created')
    return parser.parse_args()


if __name__=="__main__":
    main(parse_args())
