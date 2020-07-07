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
    parser.add_argument('--hic', type = str, help='Iteratively corrected Hi-C data')
    parser.add_argument('--chipseq', type = str, help='ChIP-Seq data')
    parser.add_argument('--np', type = int, default = 1, help='Number of processors')
    parser.add_argument('--resolution', type = int, default = 1, help='Resolution of Hi-C matrices')
    parser.add_argument('--stepsize', type = float, default = 0.5, help='Step size to increment gamma parameter')
    parser.add_argument('--gamma_max', type = float, default = 4, help='Max gamma parameter')
    parser.add_argument('--hic_format', type = str, default = 'txt.gz', help='Hi-C matrices input format for armatus')
    return parser


def check_path(path, folder_name, name):
    dirName = os.path.join(path, folder_name, name + "/")
    if not os.path.exists(dirName):
        os.makedirs(dirName, exist_ok=True)
    else:
        for old_file in glob.glob(os.path.join(dirName, "*")):
            os.remove(old_file)
    return dirName


def run_armatus(args, data, samplename):
    path = os.path.join(sys.path[0], "output")
    path_to_hic = check_path(path, "data", samplename)
    for chromosome in data.keys():
        path_to_file = os.path.join(path_to_hic, chromosome + "." + args.hic_format)
        submap = data[chromosome]
        np.savetxt(path_to_file, submap, delimiter = '\t', fmt = '%.3f')
        armatus_output = check_path(path, "tads/" + samplename, chromosome)
        caller = tadcaller.armatus(args.np, args.resolution, path_to_file, args.gamma_max, args.stepsize, armatus_output, chromosome)


def main():
    parser = get_parser()
    args = parser.parse_args()
    
    chip_data = dataloader.load_chipseq(args.chipseq)
    df = pd.DataFrame(columns = ["Gamma"])

    for filepath in glob.glob(args.hic + 'Control_Cur_0.hdf5'):
        samplename = os.path.split(filepath)[1].split(".")[0]
        print("Filename: " + samplename)
        
        hic_data = dataloader.load_hic(filepath)# load HiC map
        run_armatus(args, hic_data, samplename) # run c++ armatus

        chr_length = []
        for map in list(hic_data.values()):
            chr_length.append(map.shape[0])

        ind = tadnumeration.tadnumeration(list(hic_data.keys()), args.resolution, np.array(chr_length), samplename, args.gamma_max, args.stepsize)
        
        chip_short = chip_data[['Chr', 'Bin', samplename]]
        out = staircaller.get_stair(ind, chip_short, samplename)
        
        df_sample = pd.DataFrame(out.items(), columns = ['Gamma', samplename])
        df = pd.merge(df, df_sample, on = "Gamma", how='outer', left_index=True)
        
    visualization.plotAmplitude(df, output_path = "output/figures/StairAmplitude.png", dpi = 300)


    
if __name__ == '__main__':
    start_time = time.time()
    main()
    print('Execution time: {} sec'.format(time.time()-start_time))
