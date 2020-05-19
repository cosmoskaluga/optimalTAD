import numpy as np
import sys
import os
import glob
import h5py


def keys_parcing(obj):
    keys = list(obj.keys())
    index = keys.index('_self_key')
    keys = keys[:index]
    chromosomes = []
    n = int(np.sqrt(len(keys)))
    for i in range(n):
        chromosomes.append(keys[i*(n+1)])
    
    return chromosomes


def chr_names_parsing(obj):
    hic_chromosomes = str(obj['genomeIdxToLabel'][()])
    i = 0
    chrs = []
    while i < len(hic_chromosomes):
        if hic_chromosomes[i] == '\'':
            chr_name = ''
            i+=1
            while hic_chromosomes[i]!= '\'':
                chr_name += hic_chromosomes[i]
                i+=1
            chrs.append(chr_name)
        i+=1
    
    return chrs

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
    
    chromosomes = keys_parcing(f)
    #chromosome_names = ['2L', '2R', '3L', '3R', '4', 'X']
    chromosome_names = chr_names_parsing(f)
    matrices = {}
    
    for chrm, chrm_names in zip(chromosomes, chromosome_names):
        arr = f[chrm][()]
        
        if len(arr) > 1:
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
