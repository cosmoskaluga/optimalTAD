import lavaburst as lb
import numpy as np
import pandas as pd
import sys
import os
import glob
import h5py

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    resolution = sys.argv[3]
    
    A = np.loadtxt(input_file)
    
    for gamma in np.arange(0., 1.2, 0.1):
        S = lb.scoring.armatus_score(A, gamma = gamma)
        model = lb.SegModel(S)
    
        segments = model.optimal_segmentation()
        segments = segments[segments[:, 1] - segments[:, 0] > 2]
        segments = segments*int(resolution)
    
        chr_names = np.array(['Chr' for i in range(len(segments))])
        chr_names = pd.DataFrame(chr_names, columns = ['Chr'])
    
        segments = pd.DataFrame(segments, columns = ['Start', 'End'])
        segments[['End']] = segments[['End']] - 1
        df_row_reindex = pd.concat([chr_names, segments], ignore_index=True, axis = 1)
        
        gamma_str = str(round(gamma, 3))
        if gamma_str[2] != '0':
            gamma_str += '.0'
        path = output_file + '.gamma.' + gamma_str + '.txt'
        table_txt = df_row_reindex[::-1].to_csv(path, sep = '\t', index = None, header=False)
    

