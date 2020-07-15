import numpy as np
import pandas as pd
import h5py
import os
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
