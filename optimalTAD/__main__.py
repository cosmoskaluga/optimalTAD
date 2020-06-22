import argparse
import glob
import os
import numpy as np
from . import dataloader, tadcaller
import sys

def get_parser():
    parser = argparse.ArgumentParser(description='optimalTAD: Topologically Associating Domain optimal set prediction')
    parser.add_argument('--hic', type = str)
    parser.add_argument('--chipseq', type = str)
    parser.add_argument('--np', type = int, default = 1)
    parser.add_argument('--resolution', type = int, default = 1)
    parser.add_argument('--stepsize', default = 0.5)
    parser.add_argument('--gamma_max', default = 4)
    parser.add_argument('--hic_format', type = str, default = 'txt.gz')
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()

    for file in glob.glob(args.hic + 'Control_Cur_0.hdf5'):
        print("Filename: " + os.path.split(file)[1])
        maps = dataloader.hic(file)
        hic_data = maps.load_data()
        tadcaller.run_tadcaller(args, hic_data, file)


    
if __name__ == '__main__':
    main()
