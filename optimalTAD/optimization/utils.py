import pandas as pd
import numpy as np
import os
import glob
import shutil
import logging

log = logging.getLogger(__name__)


def optimal_gamma(data):
    samples = data.columns[1:]
    gamma = data.Gamma
    for sample in samples:
        maxindex = data[sample].idxmax()
    
    return round(gamma[maxindex], 2)


def select_optimal_tads(tads, optimal_gamma, samplename):
    dst = check_path('output/', 'optimal_gamma', samplename)
    merged_tads = pd.DataFrame(columns = ['Chr', 'Start', 'End'])
    for src in tads[optimal_gamma]:
        tad = pd.read_csv(src, header = None, names = ['Chr', 'Start', 'End'], sep = '\t')
        tad = tad[::-1]
        tad.End += 1
        merged_tads = pd.concat([merged_tads, tad])
    merged_tads.to_csv(dst + 'domains.tad', header = True, index = False)


def nan_array_comparison(func, arr, thresh):
    # https://stackoverflow.com/questions/47340000/how-to-get-rid-of-runtimewarning-invalid-value-encountered-in-greater
    # by Divakar
    bool_arr = ~np.isnan(arr)
    bool_arr[bool_arr] = func(arr[bool_arr], thresh)
    return bool_arr


def progressbar (iteration, total):
    prefix = ' '*21
    suffix = 'complete'
    length = 47
    decimals = 1
    fill = '█'
    
    percent = ('{0:.' + str(decimals) + 'f}').format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    if iteration == total:
        print()

def check_path(path, folder_name, name = None):
    dirName = os.path.join(path, folder_name, name + '/')
    if not os.path.exists(dirName):
        os.makedirs(dirName, exist_ok=True)
    else:
        for old_file in glob.glob(os.path.join(dirName, '*')):
            os.remove(old_file)
    return dirName


