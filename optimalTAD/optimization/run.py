import logging
import io
import sys
import os
import numpy as np
import pandas as pd
import cooler

from . import hicloader, chipseqloader
from . import tadcaller
from . import tadnumeration
from . import staircaller
from . import utils
from . import plotter


def main(args, cfg, log):
    df = pd.DataFrame(columns = ['Gamma'])
    stair_dict = {}
    best_gamma_array = []
    
    hic_files, chipseq_files, samplenames = utils.check_filenames(args.hic, args.chipseq)
    sample_index = 0

    args.output = utils.uniquify_path(args.output)

    for hic_path, chipseq_path, samplename in zip(hic_files, chipseq_files, samplenames):
        #samplename = os.path.split(hic_path)[1].split('.')[0]
        log.info('\033[1m' + 'Samplename: ' + samplename + '\033[0m')
    
        if sample_index < len(np.unique(hic_files)):        
            set_chromosomes = cfg.get('chromosomes', 'set_chromosomes')
            if not args.mammal:
                log.info('Load Hi-C data')
                HicLoader = hicloader.HiC(hic_path,
                                        samplename,
                                        args.hic_format,
                                        args.resolution,
                                        set_chromosomes,
                                        balance = eval(cfg.get('hic', 'balance')), 
                                        rec_size = int(cfg.get('hic', 'recommended_size_armatus'))
                                        )
                                  
                chromsize = HicLoader(args.empty_row_imputation,
                                    args.truncation,
                                    shrinkage_min = float(cfg.get('hic', 'shrinkage_min')),
                                    shrinkage_max = float(cfg.get('hic', 'shrinkage_max')),
                                    log2_hic = args.log2_hic, 
                                    output_folder = args.output)
                              
                chrs = np.array(list(chromsize.keys()), dtype = str)
                sizes = np.fromiter(chromsize.values(), dtype = int)

                log.info('Run armatus on {} chromosomes:'.format(len(chrs)))
                tadcaller.run_armatus(args, chromsize, samplename, args.output)

                log.info('Calculate indexes')
                ind, tads = tadnumeration.get_numeration(chrs, args.resolution, sizes, samplename, args.gamma_max, args.stepsize, args.output)
                blacklist = False
            else:
                log.info('Load Hi-C data and detect boundaries using IS method:')

                insulation_table = tadcaller.run_IS(path = hic_path, args = args, set_chromosomes = set_chromosomes)
                
                if args.save_insulation_score:
                    utils.save_table(insulation_table, args.output, cfg.get('output', 'path_to_insulation_score'), samplename + ".insulation_score")

                tads, blacklist = tadcaller.ins_table2tads(insulation_table, args.balance)
                ind, chromsize = tadnumeration.get_numeration_mammal(tads, args.resolution)

            sample_index += 1

        log.info('Load epigenetic data')
        ChipSeqLoader = chipseqloader.ChipSeq(chipseq_path, set_chromosomes, list(chromsize.keys()), chromsize, args.resolution)
        chip_data = ChipSeqLoader(args.log2_chip, args.zscore_chip, blacklist)
        
        if args.save_normalized_chip:
            chip_name = os.path.splitext(os.path.split(chipseq_path)[1])[0]
            utils.save_table(chip_data, args.output, cfg.get('output', 'path_to_normalized_chip'), chip_name)
            # path_to_normalized_chip = os.path.join(utils.check_path(args.output, '', cfg.get('output', 'path_to_normalized_chip')), chip_name + ".csv")
            # chip_data.to_csv(path_to_normalized_chip, header = True, index=False)
        
        log.info('Calculate amplitudes')
        stairs, amplitudes = staircaller.get_stairs(ind,
                                                    chip_data,
                                                    index_min = args.index_min,
                                                    index_max = args.index_max,
                                                    acetyl_min = cfg.getint('stair', 'acetyl_min'),
                                                    acetyl_max = cfg.getint('stair', 'acetyl_max'), 
                                                    mammals = args.mammal)
            
        #df_sample = pd.DataFrame(amplitudes.items(), columns = ['Gamma', samplename])
        df_sample = pd.DataFrame(list(amplitudes.values()), columns = ['Gamma', 
                                                         samplename, 
                                                         samplename + '_CI_lower', 
                                                         samplename + '_CI_upper'])
        if args.mammal:
            df_sample.Gamma = df_sample.Gamma.astype(int)
            df_sample = df_sample.sort_values(by=['Gamma'])
            tads = {int(k):v for k,v in tads.items()}
            stairs = {int(k):v for k,v in stairs.items()}

        best_gamma = utils.optimal_gamma(df_sample)
        log.info('The optimal gamma for {} is {}'.format(samplename, best_gamma))
        
        path_to_stairs = os.path.join(args.output, cfg.get('output', 'path_to_stair_data'))
        utils.save_stairs(stairs, 
                            index_min = args.index_min,
                            index_max = args.index_max, 
                            output_path = path_to_stairs)


        utils.select_optimal_tads(tads, best_gamma, samplename, args.output, mammal = args.mammal)

        stair_dict[samplename] = stairs[best_gamma]
        best_gamma_array.append(best_gamma)
        
        df = pd.merge(df, df_sample, on = 'Gamma', how = 'outer')
        log.info('Done!')
        print()
    
    
    path_to_amplitudes_plot =  os.path.join(args.output, cfg.get('output', 'path_to_amplitude_figure') + '.' + cfg.get('output', 'figure_postfix'))
    plotter.plotAmplitude(df, output_path = path_to_amplitudes_plot, dpi = cfg.getint('output', 'figure_dpi'))


    path_to_stairs_plot =  os.path.join(args.output, cfg.get('output', 'path_to_best_stair_figure') + '.' + cfg.get('output', 'figure_postfix'))
    stair_df = pd.DataFrame(stair_dict)
    plotter.plotStair(stair_df,
                      best_gamma_array,
                      index_min = args.index_min,
                      index_max = args.index_max,
                      output_path = path_to_stairs_plot,
                      dpi = cfg.getint('output', 'figure_dpi'), 
                      path_to_stair_dataframe = os.path.join(args.output, '', cfg.get('output', 'path_to_best_stair_data')))

    df.to_csv(os.path.join(args.output, '', cfg.get('output', 'path_to_amplitude_file')), header = True, index=False)
