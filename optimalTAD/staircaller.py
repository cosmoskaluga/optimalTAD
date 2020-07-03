import os
import sys
import numpy as np
import pandas as pd


def get_stair(index_data, chip, samplename):
    kb_list = np.arange(-4, 4, 1)
    acetyl_max = 5
    acetyl_min = -3
    dict_stairs = {key: None for key in index_data.keys()}
        
    for gamma in index_data.keys():
        chip_data = chip
        dist_data = pd.DataFrame(index_data[gamma][0], columns = ['Chr', 'Bp', 'Index'])
        chromosomes = index_data[gamma][1]
            
        interTAD = []
        TAD = []
            
        for chrm in chromosomes:
            df_dist = dist_data.loc[dist_data['Chr'] == chrm]
            df_chip = chip_data.loc[chip_data['Chr'] == chrm]
    
            for index in kb_list:
                row = df_dist.loc[df_dist['Index'] == str(index)]
                bins = row.Bp.values.astype(int)
                d_chip = df_chip.loc[df_chip['Bin'].isin(bins)]
                acetyl_val = d_chip[samplename].values
                acetyl_val = acetyl_val[~np.isnan(acetyl_val)]
                acetyl_val = acetyl_val[(acetyl_val > acetyl_min) & (acetyl_val < acetyl_max)]
                
                if index <= -1:
                    interTAD = np.append(interTAD, acetyl_val)
                if index >= 1:
                    TAD = np.append(TAD, acetyl_val)
        
            amplitude = np.median(interTAD) - np.median(TAD)
            dict_stairs[gamma] = amplitude
        
    return dict_stairs



