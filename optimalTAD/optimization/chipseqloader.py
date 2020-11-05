import numpy as np
import pandas as pd
import h5py
import os
import sys
import glob
import logging

from itertools import repeat
from . import utils


log = logging.getLogger(__name__)

bedgraph_extensions = ['.bedgraph', '.bedGraph', '.BedGraph', '.bg']
bigwig_extensions = ['.bigwig', '.bigWig', '.BigWig', '.bw']



def get_bedgraph(self):
    with open(self.path) as file_handler:
        data = np.array([])
        for line in file_handler:
            if line.startswith('chr'):
                row_content = np.array(line.split())
                if not data.size == 0:
                    data = np.vstack((data, row_content))
                else:
                    data = row_content
    
        df_chip = pd.DataFrame(data, columns = ['Chr', 'Start', 'End', 'Score'])
        df_chip = df_chip.replace('NA', 'nan')
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

    return df_chip



class ChipSeq:
    def __init__(self, path):
        self.path = path
        self.extension = os.path.splitext(path)[1]
        
        accepted_extensions = bedgraph_extensions + bigwig_extensions
        if self.extension not in accepted_extensions:
            log.error('Incompatible format of ChIP-seq file!')
            sys.exit(1)

    def __call__(self, log2_chip, set_chromosomes, zscore_chip):
        if self.extension in bedgraph_extensions:
            df_chip = get_bedgraph(self)
        else:
            df_chip = get_bigwig_file(self)
    
        if set_chromosomes != 'None':
            chrnames = set_chromosomes.split(',')
            df_chip = df_chip.loc[df_chip['Chr'].isin(chrnames)]
        
        if log2_chip:
            score = df_chip.Score.values
            bool_arr = utils.nan_array_comparison(np.less, score, 1e-10)
            score[bool_arr] = np.nan
            df_chip.Score = np.log2(df_chip.Score.values)
            df_chip = df_chip.replace(np.inf, np.nan)
        
        if zscore_chip:
            df_chip.Score = (df_chip.Score - df_chip.Score.mean()) / df_chip.Score.std(ddof=0)

        return df_chip


