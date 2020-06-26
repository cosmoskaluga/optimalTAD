import argparse
import glob
import os
import numpy as np
from . import dataloader, tadcaller, tadnumeration
import sys

def get_parser():
    parser = argparse.ArgumentParser(description='optimalTAD: Topologically Associating Domain optimal set prediction')
    parser.add_argument('--hic', type = str)
    parser.add_argument('--chipseq', type = str)
    parser.add_argument('--np', type = int, default = 1)
    parser.add_argument('--resolution', type = int, default = 1)
    parser.add_argument('--stepsize', type = float, default = 0.5)
    parser.add_argument('--gamma_max', type = float, default = 4)
    parser.add_argument('--hic_format', type = str, default = 'txt.gz')
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()

    for file in glob.glob(args.hic + 'Control_Cur_0.hdf5'):
        samplename = os.path.split(file)[1]
        print("Filename: " + samplename)
        maps = dataloader.Hic(file)
        hic_data = maps.load_data()
        tadcaller.run_tadcaller(args, hic_data, file)

        hic_maps = list(hic_data.values())
        chr_length = []
        for map in hic_maps:
            chr_length.append(map.shape[0])

        ind = tadnumeration.Indexing(list(hic_data.keys()), args.resolution, np.array(chr_length), samplename.split(".")[0], args.gamma_max, args.stepsize).get_indexes()
        


    
if __name__ == '__main__':
    main()
