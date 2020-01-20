import numpy as np
import sys
import os
import glob
import h5py

def matrix_preparation(fname, correction = True, log2transform = True, output_format = "gz"):
    """ Matrix extraction.
        
        Parameters
        ----------
        ``fname``: string
        first signal
        
        ``correction``: bool
        correction
        
        ``log2transform``: bool
        log2 transformation
        
        Returns
        -------
        class 'dict'
        HiC contact matrices
    """
    
    f = h5py.File(fname, 'r')
    
    chromosomes = ['0 0', '1 1', '2 2', '3 3', '4 4', '5 5']
    chromosome_names = ['2L', '2R', '3L', '3R', '4', 'X']
    matrices = {}
    
    for chrm, chrm_names in zip(chromosomes, chromosome_names):
        arr = f[chrm][()]
        if correction == 'True':
            arr[arr < 0.5] = 0.5
            arr[arr > 1024] = 1024
        if log2transform == 'True':
            arr = np.log2(arr)
        
        matrices[chrm_names] = arr
    
        dirName = fname[9:-5]
        if not os.path.exists(dirName):
            os.makedirs('./Output/HiCmaps/' + dirName, exist_ok=True)
        
        if output_format == "txt":
            np.savetxt('./Output/HiCmaps/' + dirName + '/' + chrm_names + '.txt', arr, delimiter = '\t', fmt = '%0.3f')
        else:
            np.savetxt('./Output/HiCmaps/' + dirName + '/' + chrm_names + '.txt.gz', arr, delimiter = '\t', fmt = '%.3f')


if __name__ == "__main__":
    correction = sys.argv[1]
    log2transform = sys.argv[2]
    
    for fname in glob.glob('HiC_data/*.hdf5'):
        print("     Extract"+" ." + fname + " file")
        matrix_preparation(fname, correction = correction, log2transform = log2transform, output_format = "gz")
