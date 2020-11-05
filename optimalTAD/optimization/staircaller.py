import os
import sys
import numpy as np
import pandas as pd


def get_stairs(index_data, chip_data, index_min = -5, index_max = 5, acetyl_min = -3, acetyl_max = 5):
    kb_list = np.arange(index_min, index_max, 1)
    gamma_range = index_data.keys()
    dict_amplitudes = {key: None for key in gamma_range}
    dict_stairs = {key: None for key in gamma_range}
    convert_chip = {'Chr': str, 'Start': int, 'End': int, 'Score': float}
    convert_dist = {'Chr': str, 'Bp': int, 'Index': str}
    
    for gamma in gamma_range:
        df_dist = pd.DataFrame(index_data[gamma][0], columns = ['Chr', 'Bp', 'Index'])
        chromosomes = index_data[gamma][1]

        
        df_dist = df_dist.loc[df_dist['Chr'].isin(chromosomes)]
        df_chip = chip_data.loc[chip_data['Chr'].isin(chromosomes)]
        
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
                
                if idx <= -1:
                    interTAD = np.append(interTAD, acetyl_val)
                if idx >= 1:
                    TAD = np.append(TAD, acetyl_val)
            else:
                median_val.append(np.nan)
        
        amplitude = np.median(interTAD) - np.median(TAD)
        dict_stairs[gamma] = np.array(median_val)
        dict_amplitudes[gamma] = amplitude
    
    return dict_stairs, dict_amplitudes



