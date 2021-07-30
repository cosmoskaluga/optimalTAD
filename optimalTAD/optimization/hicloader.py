import numpy as np
import pandas as pd
import os
import sys
import logging

from . import utils
from . import imputation

log = logging.getLogger(__name__)

accepted_extensions = ['.mcool', '.cool', '.hdf5']



def load_hdf5(path, samplename, set_chromosomes, fileformat, empty_row_imputation, truncation, shrinkage_min, shrinkage_max, log2_transformation):
    import h5py
    path = os.path.expanduser(path)
    f = h5py.File(path, 'r')
    
    if set_chromosomes == 'None':
        labels = f['chromosomeLabels'][()].astype('<U5')
    else:
        labels = set_chromosomes.split(',')
    
    path_to_output = os.path.join(sys.path[0], 'output')
    path_to_sample = utils.check_path(path_to_output, 'data', samplename)
    
    chromsize = {}
    for lb in labels:
        path_to_file = os.path.join(path_to_sample, lb + '.' + fileformat)
        matrix = f[lb][()]
        
        if empty_row_imputation:
            matrix = imputation.diagonal_interpolation(matrix)
        """
        if truncation:
            less = nan_array_comparison(np.less, matrix, shrinkage_min)
            matrix[less] = shrinkage_min
            greater = nan_array_comparison(np.greater, matrix, shrinkage_max)
            matrix[greater] = shrinkage_max
        """
        less = utils.nan_array_comparison(np.less, matrix, shrinkage_min)
        matrix[less] = shrinkage_min
        greater = utils.nan_array_comparison(np.greater, matrix, shrinkage_max)
        matrix[greater] = shrinkage_max
        
        if log2_transformation:
            matrix = np.log2(matrix)
        
        np.savetxt(path_to_file, matrix, delimiter = '\t', fmt = '%.2f')
        chromsize[lb] = matrix.shape[0]
    return chromsize


def get_coefficients(amin, amax, cmin, cmax):
    k = (cmax - cmin)/(amax - amin)
    b = cmin - k*amin
    return k, b


def load_cool(path, samplename, set_chromosomes, fileformat, balance, empty_row_imputation, truncation, shrinkage_min, shrinkage_max, log2_transformation):
    import cooler
    path_to_output = os.path.join(sys.path[0], 'output')
    path_to_sample = utils.check_path(path_to_output, 'data', samplename)
    coolfile = cooler.Cooler(path)
    
    if set_chromosomes == 'None':
        chromosomes = coolfile.chromnames
    else:
        chromosomes = set_chromosomes.split(',')
    
    chromsize = {}
    for name in chromosomes:
        path_to_file = os.path.join(path_to_sample, name + '.' + fileformat)
        matrix = coolfile.matrix(balance=balance).fetch(name)
        
        if empty_row_imputation:
            matrix = imputation.diagonal_interpolation(matrix)
        
        length = np.shape(matrix)[0]
        idx = np.arange(length)
        matrix[idx, idx] = np.repeat(np.nan, length)
        matrix[idx[:-1], idx[:-1] + 1] = np.repeat(np.nan, length - 1)
        matrix[idx[1:], idx[1:] - 1] = np.repeat(np.nan, length - 1)
        
        
        if truncation:
            #less = nan_array_comparison(np.less, matrix, shrinkage_min)
            #greater = nan_array_comparison(np.greater, matrix, shrinkage_max)
            #matrix[less] = shrinkage_min
            #matrix[greater] = shrinkage_max
            vmin = np.unique(matrix)[0]
            if vmin == 0:
                matrix[matrix == 0] = np.unique(matrix)[1]
            elif vmin > 0:
                matrix = np.clip(matrix, np.percentile(matrix, 1), np.percentile(matrix, 99))

        amin = np.nanmin(matrix)
        amax = np.nanmax(matrix)
        k,b = get_coefficients(amin, amax, 0.5, 352.7) # uncomment this
        matrix = matrix*k + b # uncomment this 
        
        if log2_transformation:
            matrix = np.log2(matrix)

        less = utils.nan_array_comparison(np.less, matrix, 0)
        if np.any(less):
            sort_data = np.unique(matrix)
            diff = min(sort_data[1:] - sort_data[:-1])
            matrix = matrix + diff # uncomment this
            #matrix = matrix - sort_data[0] + diff # comment this
    
        matrix[np.isnan(matrix, dtype=bool)] = -1 # uncomment this
        #matrix[np.isnan(matrix, dtype=bool)] = 0 # comment this

        np.savetxt(path_to_file, matrix, delimiter = '\t', fmt = '%.2f')
        chromsize[name] = matrix.shape[0]

    return chromsize



class HiC:
    def __init__(self, path, samplename, hic_format, resolution = None, set_chromosomes = None, balance = True):
        self.path = path
        self.extension = os.path.splitext(path)[1]
        self.balance = balance
        self.samplename = samplename
        self.hic_format = hic_format
        self.set_chromosomes = set_chromosomes
        
        if self.extension not in accepted_extensions:
            log.error('Incompatible format of HiC file!')
            sys.exit(1)
        
        if self.extension == '.mcool':
            suffix = '::resolutions/' + str(resolution)
            self.path += suffix

    def __call__(self, empty_row_imputation = False, truncation = False, shrinkage_min = None, shrinkage_max = None, log2_hic = False):
        if self.extension == '.hdf5':
            chromsize = load_hdf5(self.path, self.samplename, self.set_chromosomes, self.hic_format, empty_row_imputation, truncation, shrinkage_min, shrinkage_max, log2_hic)
        else:
            chromsize = load_cool(self.path, self.samplename, self.set_chromosomes, self.hic_format, self.balance, empty_row_imputation, truncation, shrinkage_min, shrinkage_max, log2_hic)

        return chromsize

