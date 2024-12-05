import numpy as np
import pandas as pd
import os
import sys
import logging

from . import utils
from . import imputation

log = logging.getLogger(__name__)

accepted_extensions = ['.mcool', '.cool', '.hdf5']


def load_hdf5(path,
            samplename,
            set_chromosomes,
            fileformat,
            rec_size,
            empty_row_imputation,
            truncation,
            shrinkage_min,
            shrinkage_max,
            log2_transformation,
            output_folder):

    """ This function uploads Hi-C maps stored in hdf5 formated files, process them (missing value imputation, truncation, log2) and writes them into .txt files

        Parameters
        ----------
        ``path`` : str
            A path to Hi-C matrix in .hdf5 format
        ``samplename`` : str
            A name of the .hdf5 file (sample name)
        ``set_chromosomes`` : str or None
            A comma-separated string of chromosome names to be included in the analysis (e.g., "chr2L,chr2R").
        ``fileformat`` : str
            An output file extension (='.txt.gz')
        ``rec_size`` : int
            Minimal of a chromosome in bins
        ``empty_row_imputation`` : bool
            If True, performs missing value imputation within empty lines of a Hi-C matrix
        ``truncation`` : bool
            If True, performs shrinkage
        ``shrinkage_min`` : int
            A min value for shrinkage (=0.5)
        ``shrinkage_max`` : int
            A max value for shrinkage (=1024)
        ``log2_transformation`` : bool
            If True, performs log2 transformation of a Hi-C matrix
        ``output_folder`` : str
            A folder name where to save the output data

        Returns
        -------
        dict
            A dictionary that stores used chromosome names together with associated sizes of Hi-C matrices in bins
    """
    import h5py
    path = os.path.expanduser(path)
    f = h5py.File(path, 'r')

    if set_chromosomes == 'None':
        labels = f['chromosomeLabels'][()].astype('<U5')
    else:
        labels_config = set_chromosomes.split(',')
        labels = utils.check_chrnames(labels_config, f['chromosomeLabels'][()].astype('<U5'))

    path_to_output = os.path.join(os.path.realpath('.'), output_folder)
    path_to_sample = utils.check_path(path_to_output, 'data', samplename)

    chromsize = {}
    for lb in labels:
        path_to_file = os.path.join(path_to_sample, lb + '.' + fileformat)
        matrix = f[lb][()]

        if np.shape(matrix)[0] > rec_size:

            if empty_row_imputation:
                matrix = imputation.diagonal_interpolation(matrix)

            if truncation:
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
    """ This function calculates coefficients for Hi-C map scaling

        Parameters
        ----------
        ``amin`` : float
            A min value of the original Hi-C map
        ``amax`` : float
            A max value of the original Hi-C map
        ``cmin`` : float
            A min value of the scaled Hi-C map
        ``cmax`` : float
            A max value of the scaled Hi-C map

        Returns
        -------
        k,b : float
            Calculated coefficients for scaling
    """
    k = (cmax - cmin)/(amax - amin)
    b = cmin - k*amin
    return k, b


def load_cool(path,
            samplename,
            set_chromosomes,
            fileformat,
            balance,
            rec_size,
            empty_row_imputation,
            truncation,
            shrinkage_min,
            shrinkage_max,
            log2_transformation,
            output_folder):

    """ This function uploads Hi-C maps stored in .cool formated files, process them (missing value imputation, truncation, log2) and writes them into .txt files

        Parameters
        ----------
        ``path`` : str
            A path to Hi-C matrix in .cool format
        ``samplename`` : str
            A name of the .cool file (sample name)
        ``set_chromosomes`` : str or None
            A comma-separated string of chromosome names to be included in the analysis (e.g., "chr2L,chr2R").
        ``balance`` : bool
            If True, uses balanced Hi-C matrices
        ``fileformat`` : str
            An output file extension (='.txt.gz')
        ``rec_size`` : int
            Minimal of a chromosome in bins
        ``empty_row_imputation`` : bool
            If True, performs missing value imputation within empty lines of a Hi-C matrix
        ``truncation`` : bool
            If True, performs shrinkage
        ``shrinkage_min`` : int
            A min value for shrinkage (=0.5)
        ``shrinkage_max`` : int
            A max value for shrinkage (=1024)
        ``log2_transformation`` : bool
            If True, performs log2 transformation of a Hi-C matrix
        ``output_folder`` : str
            A folder name where to save the output data

        Returns
        -------
        dict
            A dictionary that stores used chromosome names together with associated sizes of Hi-C matrices in bins
    """
    import cooler

    path_to_output = os.path.join(os.path.realpath('.'), output_folder)
    path_to_sample = utils.check_path(path_to_output, 'data', samplename)
    coolfile = cooler.Cooler(path)

    if set_chromosomes == 'None':
        labels = coolfile.chromnames
    else:
        labels_config = set_chromosomes.split(',')
        labels = utils.check_chrnames(labels_config, coolfile.chromnames)

    chromsize = {}
    for name in labels:
        path_to_file = os.path.join(path_to_sample, name + '.' + fileformat)
        matrix = coolfile.matrix(balance=balance).fetch(name)
        if np.shape(matrix)[0] > rec_size:

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
            k,b = get_coefficients(amin, amax, shrinkage_min, shrinkage_max) # uncomment this
            matrix = matrix*k + b # uncomment this

            if log2_transformation:
                matrix = np.log2(matrix)

            less = utils.nan_array_comparison(np.less, matrix, 0)
            if np.any(less):
                sort_data = np.unique(matrix)
                diff = min(sort_data[1:] - sort_data[:-1])
                matrix = matrix + diff # uncomment this

            matrix[np.isnan(matrix, dtype=bool)] = -1 # uncomment this

            np.savetxt(path_to_file, matrix, delimiter = '\t', fmt = '%.2f')
            chromsize[name] = matrix.shape[0]

    return chromsize



class HiC:
    def __init__(self, path, samplename, hic_format, resolution = None, set_chromosomes = None, balance = True, rec_size = 101):
        self.path = path
        self.extension = os.path.splitext(path)[1]
        self.balance = balance
        self.samplename = samplename
        self.hic_format = hic_format
        self.set_chromosomes = set_chromosomes
        self.rec_size = rec_size

        if self.extension not in accepted_extensions:
            log.error('Incompatible format of HiC file!')
            sys.exit(1)

        if self.extension == '.mcool':
            suffix = '::resolutions/' + str(resolution)
            self.path += suffix

    def __call__(self, empty_row_imputation = False, truncation = False, shrinkage_min = None, shrinkage_max = None, log2_hic = False, output_folder = 'output/'):
        if self.extension == '.hdf5':
            chromsize = load_hdf5(self.path, self.samplename, self.set_chromosomes, self.hic_format, self.rec_size, empty_row_imputation, truncation, shrinkage_min, shrinkage_max, log2_hic, output_folder)
        else:
            chromsize = load_cool(self.path, self.samplename, self.set_chromosomes, self.hic_format, self.balance, self.rec_size, empty_row_imputation, truncation, shrinkage_min, shrinkage_max, log2_hic, output_folder)

        return chromsize
