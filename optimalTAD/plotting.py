import logging
import io
import sys
import os
import numpy as np
import pandas as pd

from . visualization import hicplotter


def main(args, cfg, log):
    plotter = hicplotter.Plot(args.hic, 
    							args.region, 
    							args.resolution, 
    							args.log2_chip, 
    							args.zscore_chip, 
    							args.chipseq)

    plotter.plotHiC(cfg['hic_text'], cfg['cmap'], int(cfg['nticks']))
    plotter.plotChiPSeqTrack(cfg['chip_text'], int(cfg['fontsize']))

    plotter.plotTAD(args.tad, 
    				float(cfg['vline_linewidth']), 
    				cfg['vline_linestyle'], 
    				float(cfg['tad_linewidth']), 
    				cfg['tad_linestyle'])

    plotter.saveplot(cfg['filename'], int(cfg['dpi']))
    log.info('Done!')

