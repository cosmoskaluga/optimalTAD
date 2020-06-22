import numpy as np
import h5py
import os
import glob

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


class hic:
    def __init__(self, path):
        self.path = os.path.expanduser(path)
    
    def load_data(self):
        matrices = {}
        
        f = h5py.File(self.path, 'r')
        idx = keys_parcing(f)
        labels = chr_names_parsing(f)
        for i, lb in zip(idx, labels):
            arr = f[i][()]
            arr[arr < 0.5] = 0.5
            arr[arr > 1024] = 1024
            arr = np.log2(arr)
            matrices[lb] = arr
        return matrices


