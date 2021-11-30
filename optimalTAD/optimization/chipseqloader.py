import numpy as np
import pandas as pd
import h5py
import os
import sys
import glob
import logging

from itertools import repeat
from collections import Counter
from . import utils


log = logging.getLogger(__name__)

bedgraph_extensions = ['.bedgraph', '.bedGraph', '.BedGraph', '.bg']
bigwig_extensions = ['.bigwig', '.bigWig', '.BigWig', '.bw']



def equal_sizes(data, chromsize):
    data_size = dict(Counter(data.Chr.values))
    return data_size == chromsize



def binarize_data(data, chromsize, resolution):
    if equal_sizes(data, chromsize) != True:
        df = pd.DataFrame([])
        sizes = np.fromiter(chromsize.values(), dtype = int)
        for chrname, length in zip(np.unique(data.Chr.values), sizes):
            chr_data = data.loc[data.Chr == chrname]
            arr = []
            for i in np.arange(0, resolution*length, resolution):
                d = chr_data.loc[(chr_data.Start < i+resolution) & (chr_data.Start >= i)]
                arr.append([chrname, str(i), str(i+resolution), str(d.Score.mean(skipna=True))])

            arr = np.reshape(arr, (length, 4))
            df = pd.concat([df, pd.DataFrame(arr)])
        data = df
        data.columns = ['Chr', 'Start', 'End', 'Score']
        data.replace(['inf', '-inf'], 'nan', inplace=True)
    
    return data



def get_bedgraph(self):
    df_chip = pd.read_csv(self.path, sep = ' ', comment = 't', header = None, names = ['Chr', 'Start', 'End', 'Score'])

    if self.chrnames != ['']:
        df_chip = df_chip.loc[df_chip['Chr'].isin(self.chrnames)]

    df_chip = df_chip.replace('NA', 'nan')
    df_chip.replace(['inf', '-inf'], 'nan', inplace=True)

    df_chip = binarize_data(df_chip, self.chromsize, self.resolution)
    convert_dict = {'Chr': str, 'Start': int, 'End': int, 'Score': float}
    df_chip = df_chip.astype(convert_dict)

    return df_chip



def get_bigwig_file(self):
    import pyBigWig
    bw = pyBigWig.open(self.path)
    if bw.isBigWig() == False:
        log.error('Incompatible format of ChIP-seq file!')
        sys.exit(1)

    df_chip = pd.DataFrame([])
    for ch in bw.chroms().keys():
        intervals = bw.intervals(ch)
        df_intervals = pd.DataFrame(intervals, columns = ['Start', 'End', 'Score'])
        chr_labels = list(repeat(ch, len(df_intervals.index)))
        df_chr = pd.concat([pd.DataFrame(chr_labels, columns = ['Chr']), df_intervals], axis = 1)
        if df_chip.empty:
            df_chip = df_chr
        else:
            df_chip = pd.concat([df_chip, df_chr])

        convert_dict = {'Chr': str, 'Start': int, 'End': int, 'Score': float}
        df_chip = df_chip.astype(convert_dict)

        if self.chrnames != ['']:
            df_chip = df_chip.loc[df_chip['Chr'].isin(self.chrnames)]

        df_chip = binarize_data(df_chip, self.chromsize, self.resolution)

    return df_chip



class ChipSeq:
    def __init__(self, path, set_chromosomes, chromsize, resolution):
        self.path = path
        self.extension = os.path.splitext(path)[1]
        self.chromsize = chromsize
        self.resolution = resolution

        if set_chromosomes != 'None':
            self.chrnames = set_chromosomes.split(',')
        else:  
            self.chrnames = ['']


        accepted_extensions = bedgraph_extensions + bigwig_extensions
        if self.extension not in accepted_extensions:
            log.error('Incompatible format of ChIP-seq file!')
            sys.exit(1)

    def __call__(self, log2_chip, zscore_chip):
        if self.extension in bedgraph_extensions:
            df_chip = get_bedgraph(self)
        else:
            df_chip = get_bigwig_file(self)
        
        if log2_chip:
            score = df_chip.Score.values
            bool_arr = utils.nan_array_comparison(np.less, score, 1e-10)
            score[bool_arr] = np.nan
            df_chip.Score = np.log2(df_chip.Score.values)
            df_chip = df_chip.replace(np.inf, np.nan)
        
        if zscore_chip:
            df_chip.Score = (df_chip.Score - df_chip.Score.mean()) / df_chip.Score.std(ddof=0)

        return df_chip


