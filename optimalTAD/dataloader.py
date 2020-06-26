import numpy as np
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


