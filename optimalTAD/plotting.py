import logging
import io
import sys
import os
import numpy as np
import pandas as pd

from . visualization import hicplotter
from . optimization import utils


def main(args, cfg, log):    
    chromosome, start_bin, end_bin = utils.split_chromosome_input(args.region, args.resolution)
    path_to_hic = os.path.join('output/data', args.samplename, chromosome + '.txt.gz')

    plotter = hicplotter.Plot(path_to_hic,
                                chromosome, 
    							start_bin,
                                end_bin,
    							args.resolution, 
    							args.log2_chip, 
    							args.zscore_chip, 
    							args.chipseq)

    plotter.plotHiC(cfg['hic_text'], cfg['cmap'], int(cfg['nticks']))

    plotter.plotTAD(args.tad, float(cfg['tad_linewidth']), cfg['tad_linestyle'])

    plotter.plotChipSeqTrack(cfg['chip_text'], 
                            float(cfg['vline_linewidth']), 
                            cfg['vline_linestyle'], 
                            int(cfg['fontsize']))

    plotter.saveplot(cfg['filename'], int(cfg['dpi']))
    log.info('Done!')

