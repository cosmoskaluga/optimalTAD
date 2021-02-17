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


def split_chromosome_input(region, resolution):
    region_split = region.split(':')
    try:
        chromosome = region_split[0]
        if len(region_split) == 2:
            coordinates = region_split[1].split('-')
            start_bin = int(coordinates[0].replace(',',''))
            end_bin = int(coordinates[1].replace(',',''))

            start_bin = int(start_bin/resolution)
            end_bin = int(end_bin/resolution)
        else:
            start_bin = 0
            end_bin = None

    except ValueError:
        print('Invalid format of chromosome coordinates!')

    return chromosome, start_bin, end_bin


def get_chipname(chipnames, subfname):
    c = 0
    for name in chipnames:
        f = os.path.split(name)[1].split('.')[0]
        if f == subfname:
            break
        c+=1
    return chipnames[c]


def get_labels(start_bin, end_bin, nticks, resolution):
    interval = (end_bin - start_bin)
    ticks = list(np.linspace(0, interval, nticks).astype(int))
    pos = np.linspace(start_bin, end_bin, nticks) * resolution/1000 # kb
    pos = pos.astype(int).astype(str)
    labels = [i + ' kb' for i in pos]
    return ticks, labels


