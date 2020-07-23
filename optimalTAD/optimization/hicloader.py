import numpy as np
import pandas as pd
import os
import sys
import logging

from . import utils

log = logging.getLogger(__name__)

accepted_extensions = ['.mcool', '.cool', '.hdf5']


def load_hdf5(path, samplename, fileformat):
    import h5py
    path = os.path.expanduser(path)
    f = h5py.File(path, 'r')
    labels = f['chromosomeLabels'][()].astype('<U5')
    
    path_to_output = os.path.join(sys.path[0], 'output')
    path_to_sample = utils.check_path(path_to_output, 'data', samplename)
    
    chromsize = {}
    for lb in labels:
        path_to_file = os.path.join(path_to_sample, lb + '.' + fileformat)
        matrix = f[lb][()]
        np.savetxt(path_to_file, matrix, delimiter = '\t', fmt = '%.2f')
        chromsize[lb] = matrix.shape[0]
    
    return chromsize


def load_cool(path, samplename, fileformat, balance):
    import cooler
    path_to_output = os.path.join(sys.path[0], 'output')
    path_to_sample = utils.check_path(path_to_output, 'data', samplename)
    coolfile = cooler.Cooler(path)
    
    chromsize = {}
    for name in coolfile.chromnames:
        path_to_file = os.path.join(path_to_sample, lb + '.' + fileformat)
        chr_matrix = ccoolfile.matrix(balance=balance).fetch(name)
        np.savetxt(path_to_file, matrix, delimiter = '\t', fmt = '%.2f')
        chromsize[name] = chr_matrix.shape[0]

    return chromsize


class HiC:
    def __init__(self, path, samplename, hic_format, resolution = None, mcool_format = False, balance = False):
        self.path = path
        self.extension = os.path.splitext(path)[1]
        self.balance = balance
        self.samplename = samplename
        self.hic_format = hic_format
        
        if self.extension not in accepted_extensions:
            log.error('Incompatible format of HiC file!')
            sys.exit(1)
        
        if mcool_format == True:
            suffix = '::resolutions/' + str(resolution)
            self.path += suffix

    def __call__(self):
        if self.extension == '.hdf5':
            chromsize = load_hdf5(self.path, self.samplename, self.hic_format)
        else:
            chromsize = load_cool(self.path, self.samplename, self.hic_format, self.balance)

        return chromsize

