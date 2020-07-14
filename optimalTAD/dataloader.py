import numpy as np
import pandas as pd
import h5py
import os
import glob
import logging

log = logging.getLogger(__name__)

def load_hic(path):
    path = os.path.expanduser(path)
    matrices = {}
    f = h5py.File(path, 'r')
    labels = f['chromosomeLabels'][()].astype('<U5')
    for item in labels:
        matrices[item] = f[item][()]
    return matrices


def load_chipseq(path):
    with open(path) as file_handler:
        data = np.array([])
        for line in file_handler:
            if line.startswith("chr"):
                row_content = np.array(line.split())
                if not data.size == 0:
                    data = np.vstack((data, row_content))
                else:
                    data = row_content
                    
        df_chip = pd.DataFrame(data, columns = ['Chr', 'Start', 'End', 'Score'])
        convert_dict = {'Chr': str, 'Start': int, 'End': int, 'Score': float}
        df_chip = df_chip.astype(convert_dict)
            
    return df_chip




