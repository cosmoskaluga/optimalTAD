import numpy as np
import pandas as pd
import h5py
import os
import sys
import glob
import logging

from itertools import repeat

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

    def __call__(self):
        if self.extension in bedgraph_extensions:
            return get_bedgraph(self)
        else:
            return get_bigwig_file(self)




