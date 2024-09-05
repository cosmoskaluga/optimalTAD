import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import bootstrap
from scipy.stats import norm


def deltaH_statistic(interTAD, TAD, axis=-1):
    return np.median(interTAD, axis=axis) - np.median(TAD, axis=axis)


def get_bootstraps(data, statistic):
    rng = np.random.default_rng()
    return bootstrap(data, statistic, method='basic', random_state=rng, n_resamples=1000)


def get_stairs(index_data, df_chip, index_min = -5, index_max = 5, acetyl_min = -3, acetyl_max = 5, mammals = False):
    """ Calculating the difference in ChIP-seq signal between inter-TAD and TAD regions.
        
        Parameters
        ----------
        ``index_data`` : dataframe
            A dataframe consisting of the distances (=indexes) to the nearest TAD boundaries from each bin.
        ``df_chip`` : dataframe
            A dataframe with corresponding ChIP-seq signal 
        ``index_min`` : int
            Minimal distance value to consider while calculating median of ChIP-seq signal
        ``index_max`` : int
            Maximal distance value to consider while calculating median of ChIP-seq signal
        ``acetyl_min`` : int
            A lower threshold for a scaled ChIP-seq signal 
        ``acetyl_max`` : int
            An upper threshold for a scaled ChIP-seq signal 
        ``mammals`` : bool
            Telling the algorithm that the data corresponds to mammalian chromatin
        
        Returns
        -------
        ``dict_stairs`` : dict
            Median ChIP-seq values per each index used (per each distance to the nearest TAD boundary)
        ``dict_amplitudes`` : dict
            Difference in median ChIP-seq value between inter-TADs and TADs for each value of the optimized parameter 
    """

    kb_list = np.arange(index_min, index_max, 1)
    gamma_range = index_data.keys()
    dict_amplitudes = {key: None for key in gamma_range}
    dict_stairs = {key: None for key in gamma_range}
    convert_chip = {'Chr': str, 'Start': int, 'End': int, 'Score': float}
    convert_dist = {'Chr': str, 'Bp': int, 'Index': str}
    
    for gamma in gamma_range:
        df_dist = pd.DataFrame(index_data[gamma][0], columns = ['Chr', 'Bp', 'Index'])
        chromosomes = index_data[gamma][1]

        #df_dist = df_dist.loc[df_dist['Chr'].isin(chromosomes)]
        #df_chip = chip_data.loc[chip_data['Chr'].isin(chromosomes)]
        
        df_dist = df_dist.astype(convert_dist)
        df_chip = df_chip.astype(convert_chip)
        
        df_dist = df_dist.sort_values(by = ["Chr", "Bp"])
        df_chip = df_chip.sort_values(by = ["Chr", "Start"])
        
        df_dist.index = np.arange(df_dist.shape[0])
        df_chip.index = np.arange(df_chip.shape[0])
        
        median_val = []
        interTAD = []
        TAD = []

        for idx in kb_list:
            row = df_dist.loc[df_dist.Index == str(idx)]
            if row.empty == False:
                bins = row.index.values
                d_chip = df_chip.loc[df_chip.index.isin(bins)]
                acetyl_val = d_chip['Score'].values
                acetyl_val = acetyl_val[~np.isnan(acetyl_val)]
                acetyl_val = acetyl_val[(acetyl_val > acetyl_min) & (acetyl_val < acetyl_max)]
                
                if len(acetyl_val) == 0:
                    median_val.append(np.nan)
                else:
                    median_val.append(np.median(acetyl_val))

                if mammals:
                    border_bin = 1
                else: 
                    border_bin = 0

                if idx < border_bin:
                    interTAD = np.append(interTAD, acetyl_val)
                if idx >= border_bin:
                    TAD = np.append(TAD, acetyl_val)
            else:
                median_val.append(np.nan)
        
        amplitude = deltaH_statistic(interTAD, TAD)
        res = get_bootstraps((interTAD, TAD), deltaH_statistic)
        dict_stairs[gamma] = np.array(median_val)
        dict_amplitudes[gamma] = [gamma, amplitude, res.confidence_interval.low, res.confidence_interval.high]
    
    return dict_stairs, dict_amplitudes

