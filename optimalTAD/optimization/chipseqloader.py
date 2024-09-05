import numpy as np
import pandas as pd
import h5py
import os
import sys
import glob
import logging
import csv

from itertools import repeat
from collections import Counter, OrderedDict
import collections
from . import utils


log = logging.getLogger(__name__)

bedgraph_extensions = ['.bedgraph', '.bedGraph', '.BedGraph', '.bg']
bigwig_extensions = ['.bigwig', '.bigWig', '.BigWig', '.bw']



def add_chr(dic):
    """ Adding 'chr' prefix to chromosome names.
        
        Parameters
        ----------
        ``dic`` : dictionary
            A dictionary containing chromosome-size pairs
        
        Returns
        -------
        dict
            A dictionary with corrected chromosome names (chr1, chr2, etc)
    """
    corrected_names = ['chr'+str(s) for s in dic.keys()]
    return {k: v for k, v in zip(corrected_names, list(dic.values()))}



def equal_sizes(data, chromsize):
    """ Checking whether bin sizes are equal in Hi-C and ChIP-seq data.
        
        Parameters
        ----------
        ``data`` : dataframe
            ChIP-seq signal
        ``chromsize`` : dictionary
            Chromosome sizes obtained from Hi-C data
        
        Returns
        -------
        bool
            True or False
    """
    data_size = dict(Counter(data.Chr.values))
    bool_data = np.all(utils.check_prefix(data_size.keys()))
    bool_chromsize = np.all(utils.check_prefix(chromsize.keys()))
    if bool_data == bool_chromsize:
        pass
    elif not bool_data:
        data_size = add_chr(data_size)
    elif not bool_chromsize:
        chromsize = add_chr(chromsize)

    data_size = OrderedDict(sorted(data_size.items()))
    chromsize = OrderedDict(sorted(chromsize.items()))

    return data_size == chromsize



def bin_calculation(chrdata, length, chrname, resolution):
    """ ChIP-seq signal re-aggregation within a given genomic interval (=Hi-C map resolution).
        
        Parameters
        ----------
        ``chrdata`` : dataframe
            ChIP-seq signal
        ``length`` : int
            Chromosome length
        ``chrname`` : str
            Chromosome name
        ``resolution`` : int
            Desired bin size (equal to Hi-C map resolution)
        
        Returns
        -------
        np.ndarray
            Binarized ChIP-seq signal for a chromosome 
    """
    df = pd.DataFrame([])
    arr = []
    i = 0
    start_value = 0
    while i < length*resolution:
        d = chrdata.loc[(chrdata.Start < i+resolution) & (chrdata.Start >= i)]
        end_values = d.End.values
        if len(end_values) > 0:
            last_interval = end_values[-1] - (i+resolution)
            nbins = (end_values[-1] - d.Start.values[0])//resolution
            if nbins > 1:
                last_interval = last_interval%resolution
                new_bin = d.Score.values[-1]
            else:
                end_values[-1] -= last_interval
                interval_len = end_values - d.Start.values
                binvals = d.Score.values*interval_len/resolution
                new_bin = np.nansum(np.append(binvals, start_value))
                nbins = 1
            start_value = d.Score.values[-1]*last_interval/resolution
        
        else:
            new_bin = np.nan
            if i < chrdata.Start.values[0]:
                nbins = (chrdata.Start.values[0] - i)//resolution
            else:
                nbins = (length*resolution - i)//resolution
        
        for bin_idx in range(nbins):
            arr.append([chrname, str(i+bin_idx*resolution), str(i+resolution*(1+bin_idx)), str(new_bin)])
        
        i += resolution*nbins
    
    arr = np.reshape(arr, (i//resolution, 4))
    
    return arr



def binarize_data(data, chromsize, resolution):
    """ ChIP-seq signal re-aggregation within a given genomic interval (=Hi-C map resolution).
        
        Parameters
        ----------
        ``data`` : dataframe
            ChIP-seq signal
        ``chromsize`` : dictionary
            Chromosome sizes obtained from Hi-C data
        ``resolution`` : int
            Desired bin size (equal to Hi-C map resolution)
        
        Returns
        -------
        dataframe
            Binarized ChIP-seq signal
    """
    if equal_sizes(data, chromsize) != True:
        df = pd.DataFrame([])
        sizes = np.fromiter(chromsize.values(), dtype = int)
        for chrname, length in zip(np.unique(data.Chr.values), sizes):
            chr_data = data.loc[data.Chr == chrname]
            #arr = []
            #for i in np.arange(0, resolution*length, resolution):
            #    d = chr_data.loc[(chr_data.Start < i+resolution) & (chr_data.Start >= i)]
            #    arr.append([chrname, str(i), str(i+resolution), str(d.Score.mean(skipna=True))])
            
            #arr = np.reshape(arr, (length, 4))
            arr = bin_calculation(chr_data, length, chrname, resolution)
            df = pd.concat([df, pd.DataFrame(arr)])
        data = df
        data.columns = ['Chr', 'Start', 'End', 'Score']
        data.replace(['inf', '-inf'], 'nan', inplace=True)
    
    return data



def blacklist_subtraction(signal_data, bklst):
    """ Assigning nan values to blacklist bins.
        
        Parameters
        ----------
        ``signal_data`` : dataframe
            ChIP-seq signal
        ``blklst`` : dataframe
            Dataframe consisting of blacklist regions coordinates 
        
        Returns
        -------
        dataframe
            ChIP-seq signal with nans for blacklist regions
    """

    keys = list(bklst.columns.values)
    i1 = signal_data.set_index(keys).index
    i2 = bklst.set_index(keys).index
    signal_data.loc[i1.isin(i2), 'Score'] = np.nan
    return signal_data



def get_bedgraph(self, blacklist_regions = False):
    """ Read data from bedgraph file 
        
        Parameters
        ----------
        ``blacklist_regions`` : dataframe
            Dataframe consisting of blacklist regions coordinates (default False)
        
        Returns
        -------
        dataframe
            Processed and binarized (if needed) ChIP-seq signal
    """
    with open(self.path) as f:
        header = csv.Sniffer().has_header(f.read(1024))

    df_chip = pd.read_csv(self.path, sep = '\s+', comment = 't', header = None, skiprows=header*1, names = ['Chr', 'Start', 'End', 'Score'],  
                          dtype = {"Chr": object, "Start": np.int64, "End": np.int64, "Score": np.float64})  

    labels = utils.check_chrnames(self.chrnames, np.unique(df_chip.Chr))

    df_chip = df_chip.loc[df_chip['Chr'].isin(labels)]

    df_chip = df_chip.replace('NA', 'nan')
    df_chip.replace(['inf', '-inf'], 'nan', inplace=True)

    df_chip = binarize_data(df_chip, self.chromsize, self.resolution)
    convert_dict = {'Chr': str, 'Start': int, 'End': int, 'Score': float}
    df_chip = df_chip.astype(convert_dict)

    if isinstance(blacklist_regions, pd.DataFrame):
        df_chip = blacklist_subtraction(df_chip, blacklist_regions)

    return df_chip



def get_bigwig_file(self, blacklist_regions = False):
    """ Read data from bigwig file 
        
        Parameters
        ----------
        ``blacklist_regions`` : dataframe
            Dataframe consisting of blacklist regions coordinates (default False)
        
        Returns
        -------
        dataframe
            Processed and binarized (if needed) ChIP-seq signal
    """    
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
        

    labels = utils.check_chrnames(self.chrnames, np.unique(df_chip.Chr))
    df_chip = df_chip.loc[df_chip['Chr'].isin(labels)]
    df_chip = binarize_data(df_chip, self.chromsize, self.resolution)


    if isinstance(blacklist_regions, pd.DataFrame):
        df_chip = blacklist_subtraction(df_chip, blacklist_regions)

    return df_chip



class ChipSeq:
    def __init__(self, path, set_chromosomes, hic_chromosomes, chromsize, resolution):
        self.path = path
        self.extension = os.path.splitext(path)[1]
        self.chromsize = chromsize
        self.resolution = resolution

        if set_chromosomes != 'None':
            self.chrnames = set_chromosomes.split(',')
        else:  
            self.chrnames = hic_chromosomes

        accepted_extensions = bedgraph_extensions + bigwig_extensions
        if self.extension not in accepted_extensions:
            log.error('Incompatible format of ChIP-seq file!')
            sys.exit(1)

    def __call__(self, log2_chip, zscore_chip, blacklist_regions = False):
        if self.extension in bedgraph_extensions:
            df_chip = get_bedgraph(self, blacklist_regions)
        else:
            df_chip = get_bigwig_file(self, blacklist_regions)

        score = df_chip.Score.values.astype(float)
        bool_arr = utils.nan_array_comparison(np.less, score, 1e-10)
        score[bool_arr] = np.nan
        df_chip.Score = score
        
        if log2_chip:
            score = df_chip.Score.values
            df_chip.Score = np.log2(score)
            df_chip = df_chip.replace(np.inf, np.nan)
        
        if zscore_chip:
            df_chip.Score = (df_chip.Score - df_chip.Score.mean()) / df_chip.Score.std(ddof=0)

        return df_chip


