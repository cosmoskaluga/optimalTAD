import logging
import io
import sys
import os
import numpy as np
import pandas as pd

from . visualization import hicplotter


def main(args, log):
    plotter = hicplotter.Plot(args.hic, args.tad, args.region, args.resolution, args.chipseq)
    plotter.plotHiC()
    plotter.plotChiPSeqTrack()
    plotter.plotTAD()
    plotter.saveplot('output/figures/TADandChIPseq.png', dpi = 300)

