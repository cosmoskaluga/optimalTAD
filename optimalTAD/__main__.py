import argparse
import logging
import glob
from tqdm import tqdm
import io
import sys
import os
import time
import numpy as np
import pandas as pd

from . import dataloader
from . import logger
from . import tadcaller
from . import tadnumeration
from . import staircaller
from . import utils
from . import visualization

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


def run_armatus(args, data, samplename, log):
    path = os.path.join(sys.path[0], "output")
    path_to_hic = utils.check_path(path, "data", samplename)
    
    log.info("Running armatus in serial on all {} chromosomes:".format(len(data.keys())))
    num = 0
    for chromosome in data.keys():
        utils.progressbar(num, len(data.keys()))
        path_to_file = os.path.join(path_to_hic, chromosome + "." + args.hic_format)
        submap = data[chromosome]
        np.savetxt(path_to_file, submap, delimiter = '\t', fmt = '%.2f')
        armatus_output = utils.check_path(path, "tads/" + samplename, chromosome)
        caller = tadcaller.armatus(args.np, args.resolution, path_to_file, args.gamma_max, args.stepsize, armatus_output, chromosome)
        num+=1
    utils.progressbar(num, len(data.keys()))


def main():
    parser = get_parser()
    args = parser.parse_args()
    
    log = logger.initialize_logger()
    
    chip_data = dataloader.load_chipseq(args.chipseq)
    df = pd.DataFrame(columns = ["Gamma"])
    log.info("Loading ChIP-Seq data")

    for filepath in glob.glob(args.hic + '*.hdf5'):
        samplename = os.path.split(filepath)[1].split(".")[0]
        print('')
        log.info('\033[1m' + "Filename: " + samplename + '\033[0m')
        
        log.info("Loading Hi-C data")
        hic_data = dataloader.load_hic(filepath)# load HiC map
        run_armatus(args, hic_data, samplename, log) # run c++ armatus

        chr_length = []
        for map in list(hic_data.values()):
            chr_length.append(map.shape[0])
        
        log.info("Calculating indexes")
        ind = tadnumeration.get_numeration(list(hic_data.keys()), args.resolution, np.array(chr_length), samplename, args.gamma_max, args.stepsize)
        
        log.info("Calculating stair amplitude")
        chip_short = chip_data[['Chr', 'Bin', samplename]]
        out = staircaller.get_stairs(ind, chip_short, samplename)
        
        df_sample = pd.DataFrame(out.items(), columns = ['Gamma', samplename])
        df = pd.merge(df, df_sample, on = "Gamma", how='outer', left_index=True)
        log.info("Done!")

    print('')
    utils.optimal_gamma(df)
    visualization.plotAmplitude(df, output_path = "output/figures/StairAmplitude.png", dpi = 300)


    
if __name__ == '__main__':
    start_time = time.time()
    main()
    print('Execution time: {} sec'.format(time.time()-start_time))
