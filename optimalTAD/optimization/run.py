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
from . import plotter


def run_armatus(args, chromsize, samplename):
    path = os.path.join(sys.path[0], 'output')
    path_to_hic = os.path.join(path, 'data', samplename + '/')
    num = 0
    for chromosome in chromsize.keys():
        utils.progressbar(num, len(chromsize.keys()))
        path_to_file = os.path.join(path_to_hic, chromosome + '.' + args.hic_format)
        armatus_output = utils.check_path(path, 'tads/' + samplename, chromosome)
        caller = tadcaller.armatus(args.np,
                                   args.resolution,
                                   path_to_file,
                                   args.gamma_max,
                                   args.stepsize,
                                   armatus_output,
                                   chromosome)
        num += 1
    utils.progressbar(num, len(chromsize.keys()))


def main(args, cfg, log):
    df = pd.DataFrame(columns = ['Gamma'])
    stair_dict = {}
    best_gamma_array = []
    
    for hic_path, chipseq_path in zip(args.hic, args.chipseq):
        samplename = os.path.split(hic_path)[1].split('.')[0]
        log.info('\033[1m' + 'Samplename: ' + samplename + '\033[0m')
            
        log.info('Loading Hi-C data')
        set_chromosomes = cfg.get('chromosomes', 'set_chromosomes')
        HicLoader = hicloader.HiC(hic_path,
                                  samplename,
                                  args.hic_format,
                                  args.resolution,
                                  set_chromosomes,
                                  balance = eval(cfg.get('hic', 'balance')))
                                  
        chromsize = HicLoader(args.empty_row_imputation,
                              args.truncation,
                              shrinkage_min = float(cfg.get('hic', 'shrinkage_min')),
                              shrinkage_max = float(cfg.get('hic', 'shrinkage_max')),
                              log2_hic = args.log2_hic)
                              
        chrs = np.array(list(chromsize.keys()), dtype = str)
        sizes = np.fromiter(chromsize.values(), dtype = int)
        
        
        log.info('Loading epigenetic data')
        ChipSeqLoader = chipseqloader.ChipSeq(chipseq_path)
        chip_data = ChipSeqLoader(args.log2_chip, set_chromosomes, args.zscore_chip)
        
        
        log.info('Running armatus in serial on all {} chromosomes:'.format(len(chrs)))
        run_armatus(args, chromsize, samplename)
        
        
        log.info('Calculating indexes')
        ind, tads = tadnumeration.get_numeration(chrs, args.resolution, sizes, samplename, args.gamma_max, args.stepsize)
        
        #print(tads[0.0])
        
        log.info('Calculating stair amplitude')
        stairs, amplitudes = staircaller.get_stairs(ind,
                                                    chip_data,
                                                    index_min = cfg.getint('stair', 'index_min'),
                                                    index_max = cfg.getint('stair', 'index_max'),
                                                    acetyl_min = cfg.getint('stair', 'acetyl_min'),
                                                    acetyl_max = cfg.getint('stair', 'acetyl_max'))
            
        df_sample = pd.DataFrame(amplitudes.items(), columns = ['Gamma', samplename])
        best_gamma = utils.optimal_gamma(df_sample)
        log.info('The optimal gamma for {} is {}'.format(samplename, best_gamma))
        
        utils.select_optimal_tads(tads, best_gamma, samplename)
        
        stair_dict[samplename] = stairs[best_gamma]
        best_gamma_array.append(best_gamma)
        
        df = pd.merge(df, df_sample, on = 'Gamma', how = 'outer')
        log.info('Done!')
        print()
    
    
    path = cfg.get('output', 'path_to_amplitude_figure') + '.' + cfg.get('output', 'figure_postfix')
    plotter.plotAmplitude(df, output_path = path, dpi = cfg.getint('output', 'figure_dpi'))


    stair_df = pd.DataFrame(stair_dict)
    plotter.plotStair(stair_df,
                      best_gamma_array,
                      index_min = cfg.getint('stair', 'index_min'),
                      index_max = cfg.getint('stair', 'index_max'),
                      output_path = cfg.get('output', 'path_to_stair_figure') + '.' + cfg.get('output', 'figure_postfix'),
                      dpi = cfg.getint('output', 'figure_dpi'), 
                      path_to_stair_dataframe = cfg.get('output', 'path_to_stair_dataframe'))

    df.to_csv(cfg.get('output', 'path_to_amplitude_file'), header = True, index=False)
