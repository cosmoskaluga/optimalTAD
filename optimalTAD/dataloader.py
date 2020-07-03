import numpy as np
import pandas as pd
import h5py
import os
import glob


class Hic:
    def __init__(self, path):
        self.path = os.path.expanduser(path)
    
    def load_data(self):
        matrices = {}
        f = h5py.File(self.path, 'r')
        labels = f['chromosomeLabels'][()].astype('<U2')
        for item in labels:
            matrices[item] = f[item][()]
        return matrices


class ChIP_Seq:
    def __init__(self, path):
        self.path = path
    def get_data(self):
        chip_data = pd.read_csv(self.path, sep = '\t')
        return chip_data


