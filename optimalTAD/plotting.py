import logging
import io
import sys
import os
import numpy as np
import pandas as pd

from . visualization import hicplotter


def main(args, log):
    plotter = hicplotter.Plot(args.hic, args.region, args.resolution, args.log2_chip, args.zscore_chip, args.chipseq)
    plotter.plotHiC(text = 'Hi-C', cmap = 'coolwarm', nticks = 5)
    plotter.plotChiPSeqTrack(text = 'RNA-seq', fontsize = 12)
    plotter.plotTAD(args.tad, vline_linewidth = 1., vline_linestyle = 'dashed', tad_linewidth = 1.2, tad_linestyle = '--')
    plotter.saveplot('output/figures/TADandChIPseq.png', dpi = 300)
    log.info('Done!')

