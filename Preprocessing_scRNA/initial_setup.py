import os
import sys
import argparse


def main(args):
    # first check if the folder already exists
    path = os.path.join(args.work_dir, args.id, "readme.md")
    print(path)
    if os.path.exists(path):
        sys.exit("The directory already exists.")
    else:
        dir_path = os.path.join(args.work_dir, args.id)
        # os.mkdir(dir_path)
        readme = open(os.path.join(dir_path, "readme.md"), "w")
        print("#Documentation", file=readme)
        readme.close()

def parse_args():
    parser = argparse.ArgumentParser(description='Initial setup for processing new scRNAseq data')
    parser.add_argument('--id', required=True, help='Input id as a number. For example: 1, 2, 3, 4, etc...')
    parser.add_argument('--work_dir', required=True, help='Input working directory where the folder is going to be created')
    return parser.parse_args()


if __name__=="__main__":
    main(parse_args())
