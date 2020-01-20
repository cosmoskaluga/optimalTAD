import numpy as np
import pandas as pd
import sys
import os
import glob

import warnings
warnings.filterwarnings("ignore")

def numerating(l, flag):
    if flag == 'intertad':
        if l%2 == 0:
            m = list(range(-1, -1-int(l/2), -1))
            m = m + m[::-1]
        else:
            m = list(range(-1, -1-int((l-1)/2), -1))
            m = m + [-1-int((l-1)/2)] + m[::-1]
    else:
        if l%2 == 0:
            m = list(range(0, int(l/2)))
            m = m + m[::-1]
        else:
            m = list(range(0, int((l-1)/2)))
            m = m + [int(l/2)] + m[::-1]

    return(np.array(m))

def distance_calculating(d, resolution, size):
    first_col = np.array([])
    second_col = np.array([])
    length = int(size*resolution)
    bp = 0
    for row in d:
        start = row[1]
        end = row[2] + 1
        diff = end - start
        if bp != start:
            l = (start - bp)/resolution
            col1 = np.arange(bp, start - 0.01, resolution)
            col2 = numerating(l, 'intertad')
            first_col = np.append(first_col, col1)
            second_col = np.append(second_col, col2)
        
        
        l = (end - start)/resolution
        col1 = np.arange(start, end, resolution)
        col2 = numerating(l, 'tad')
        first_col = np.append(first_col, col1)
        second_col = np.append(second_col, col2)
        
        bp = int(end)

    if (length > bp):
        start = bp
        end = length
        diff = end - start
        l = (end - bp)/resolution
        col1 = np.arange(bp, end - 0.01, resolution)
        col2 = numerating(l, 'intertad')
        first_col = np.append(first_col, col1)
        second_col = np.append(second_col, col2)

    #chrm_list = ['chr2L' for i in range(size)]
    out = np.vstack((chrm_list, first_col, second_col))
    return np.transpose(out)


if __name__ == "__main__":
    dirName = sys.argv[1][12:]
    
    chrm = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chrX']
    chrm_index = [0, 1151, 2209, 3437, 4833, 6023 - 68]

    chrm_list = []
    for ind in range(len(chrm)):
        a = [chrm[ind] for i in range(chrm_index[ind+1] - chrm_index[ind])]
        chrm_list+=a

    resolution = 20000

    chromosomeStarts = np.array([0, 1151, 2209, 3437, 4833])
    chromosomeStarts *= resolution

    if not os.path.exists(dirName):
        os.makedirs('./Output/NumeratedBorders/' + dirName, exist_ok=True)
    
    for f2L, f2R, f3L, f3R, fX in zip(glob.glob('Output/TADs/' + dirName + '/2L*.txt'),
                                      glob.glob('Output/TADs/' + dirName + '/2R*.txt'),
                                      glob.glob('Output/TADs/' + dirName + '/3L*.txt'),
                                      glob.glob('Output/TADs/' + dirName + '/3R*.txt'),
                                      glob.glob('Output/TADs/' + dirName + '/X*.txt')):
    
        data_2L = pd.read_csv(f2L, header = None, names = ['Chr', 'start', 'end'], sep = '\t')
        data_2R = pd.read_csv(f2R, header = None, names = ['Chr', 'start', 'end'], sep = '\t')
        data_3L = pd.read_csv(f3L, header = None, names = ['Chr', 'start', 'end'], sep = '\t')
        data_3R = pd.read_csv(f3R, header = None, names = ['Chr', 'start', 'end'], sep = '\t')
        data_X = pd.read_csv(fX, header = None, names = ['Chr', 'start', 'end'], sep = '\t')
    
        data_2L.Chr = pd.DataFrame(['Chr2L' for i in range(len(data_2L))])
        data_2R.Chr = pd.DataFrame(['Chr2R' for i in range(len(data_2R))])
        data_3L.Chr = pd.DataFrame(['Chr3L' for i in range(len(data_3L))])
        data_3R.Chr = pd.DataFrame(['Chr3R' for i in range(len(data_3R))])
        data_X.Chr = pd.DataFrame(['ChrX' for i in range(len(data_X))])

        for arr, i in zip([data_2L, data_2R, data_3L, data_3R, data_X], range(5)):
            arr.start = arr.start + chromosomeStarts[i]
            arr.end = arr.end + chromosomeStarts[i]

        data = [data_2L, data_2R, data_3L, data_3R, data_X]
        data_conc = pd.concat([data_2L[::-1], data_2R[::-1], data_3L[::-1], data_3R[::-1], data_X[::-1]])
        data_conc.index = range(0, len(data_conc))

        array = distance_calculating(data_conc.values, resolution, 6023-68)
        df = pd.DataFrame(array, columns = ['Chr', 'Bp', 'Number'])

        export_csv = df.to_csv('./Output/NumeratedBorders/'+f2L[12:-4]+'.dist', index = None, header=False)







