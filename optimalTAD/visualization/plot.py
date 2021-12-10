import logging
import io
import sys
import os
import numpy as np
import pandas as pd

from . import settings, dataloader
from .. optimization import utils

log = logging.getLogger(__name__)


def main(args, cfg, log):    
    chromosome, start_bin, end_bin = utils.split_chromosome_input(args.region, args.resolution)
    hic_filename = os.path.join('output/data', args.samplename, chromosome + '.txt.gz')


    plotter = settings.Plot(hic_filename,
                                chromosome, 
    							start_bin,
                                end_bin,
    							args.resolution)

    plotter.plotHiC(cfg['hic_text'], cfg['cmap'], int(cfg['nticks']))

    tads = dataloader.get_domains(args.samplename, chromosome)
    plotter.plotTAD(tads, float(cfg['tad_linewidth']), cfg['tad_linestyle'])

    if args.chipseq != False:
        chromsize = plotter.get_chromsize()

        chipseq = dataloader.get_chipseq(args.chipseq, chromosome, chromsize, args.resolution, args.log2_chip, args.zscore_chip)
        plotter.plotTrack(chipseq, cfg['chip_text'], 
                                    float(cfg['vline_linewidth']), 
                                    cfg['vline_linestyle'], 
                                    int(cfg['fontsize']))

    if args.rnaseq != 'False':
        rnaseq = dataloader.get_rnaseq(args.rnaseq, chromosome, start_bin, end_bin)
        plotter.plotTrack(rnaseq, cfg['rnaseq_text'], 
                                float(cfg['vline_linewidth']), 
                                cfg['vline_linestyle'], 
                                int(cfg['fontsize']))


    plotter.show()
    plotter.saveplot(cfg['filename'], int(cfg['dpi']))
    log.info('Done!')

