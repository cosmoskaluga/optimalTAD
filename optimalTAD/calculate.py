import logging
import io
import sys
import os
import numpy as np
import pandas as pd

from . import hicloader, chipseqloader
from . import tadcaller
from . import tadnumeration
from . import staircaller
from . import utils
from . import visualization


def run_armatus(args, chromsize, samplename):
    path = os.path.join(sys.path[0], 'output')
    path_to_hic = os.path.join(path, 'data', samplename + '/')
    num = 0
    for chromosome in chromsize.keys():
        utils.progressbar(num, len(chromsize.keys()))
        path_to_file = os.path.join(path_to_hic, chromosome + '.' + args.hic_format)
        armatus_output = utils.check_path(path, 'tads/' + samplename, chromosome)
        caller = tadcaller.armatus(args.np, args.resolution, path_to_file, args.gamma_max, args.stepsize, armatus_output, chromosome)
        num+=1
    utils.progressbar(num, len(chromsize.keys()))


def main(args, log):
    df = pd.DataFrame(columns = ['Gamma'])
    
    for hic_path, chipseq_path in zip(args.hic, args.chipseq):
        samplename = os.path.split(hic_path)[1].split('.')[0]
        log.info('\033[1m' + 'Samplename: ' + samplename + '\033[0m')
            
        log.info('Loading Hi-C data')
        HicLoader = hicloader.HiC(hic_path, samplename, args.hic_format, mcool_format = False)
        chromsize = HicLoader()
        chrs = np.array(list(chromsize.keys()), dtype = str)
        sizes = np.fromiter(chromsize.values(), dtype = int)
            
        log.info('Loading ChIP-seq data')
        ChipSeqLoader = chipseqloader.ChipSeq(chipseq_path)
        chip_data = ChipSeqLoader()
            
        log.info('Running armatus in serial on all {} chromosomes:'.format(len(chrs)))
        run_armatus(args, chromsize, samplename)
            
        log.info('Calculating indexes')
        ind, tads = tadnumeration.get_numeration(chrs, args.resolution, sizes, samplename, args.gamma_max, args.stepsize)
            
        log.info('Calculating stair amplitude')
        out = staircaller.get_stairs(ind, chip_data)
            
        df_sample = pd.DataFrame(out.items(), columns = ['Gamma', samplename])
        gamma_best = utils.optimal_gamma(df_sample)
        log.info('The optimal gamma for {} is {}'.format(samplename, gamma_best))
            
        df = pd.merge(df, df_sample, on = 'Gamma', how='outer', left_index=True)
        log.info('Done!')
        print('')
        
    visualization.plotAmplitude(df, output_path = 'output/figures/StairAmplitude.png', dpi = 300)
    df.to_csv('output/amplitudes.csv', header = True, index=False)
