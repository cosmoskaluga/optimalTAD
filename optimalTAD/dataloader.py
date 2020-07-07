import numpy as np
import pandas as pd
import h5py
import os
import glob

def load_hic(path):
    path = os.path.expanduser(path)
    matrices = {}
    f = h5py.File(path, 'r')
    labels = f['chromosomeLabels'][()].astype('<U2')
    for item in labels:
        matrices[item] = f[item][()]
    return matrices

def load_chipseq(path):
    chip_data = pd.read_csv(path, sep = '\t')
    return chip_data




