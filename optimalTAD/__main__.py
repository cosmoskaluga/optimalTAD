import argparse
import glob
import os
import time
import numpy as np
import pandas as pd
from . import dataloader, tadcaller, tadnumeration
from . import staircaller
from . import visualization
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
    
    chip_data = dataloader.ChIP_Seq(args.chipseq).get_data()

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

        sample = samplename.split(".")[0]
        ind = tadnumeration.Indexing(list(hic_data.keys()), args.resolution, np.array(chr_length), sample, args.gamma_max, args.stepsize).get_indexes()
        
        ship_short = chip_data[['Chr', 'Bin', sample]]
        out = staircaller.get_stair(ind, ship_short, sample)
        visualization.plotAmplitude(out, sample, output_path = "output/figures/StairAmplitude.png", dpi = 300)



    
if __name__ == '__main__':
    start_time = time.time()
    main()
    print('Execution time: {} sec'.format(end_time-time.time()))
