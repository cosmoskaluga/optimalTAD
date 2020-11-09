import numpy as np
import pandas as pd
import os
import logging

from .. optimization import chipseqloader

log = logging.getLogger(__name__)

def get_domains(samplename, chromosome):
    path = os.path.join('output/optimal_gamma', samplename, 'domains.tad')
    df = pd.read_csv(path)
    df = df.loc[df.Chr == chromosome]
    return df

def get_chipseq(samplename, chromosome, log2_chip = False, zscore_chip = False):
	ChipSeqLoader = chipseqloader.ChipSeq(samplename)
	chip_data = ChipSeqLoader(log2_chip, chromosome, zscore_chip)
	chip_data = chip_data.loc[chip_data.Chr.isin([chromosome])]
	return chip_data

def get_rnaseq(samplename, chromosome, start_bin, end_bin):
	df = pd.read_csv(samplename, sep = ' ', header = None, names = ['Chr', 'Start', 'End', 'Score'])
	df = df.loc[df.Chr == chromosome]
	return df