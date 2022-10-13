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


def save_stairs(data, index_min, index_max, output_path):
    data = pd.DataFrame(data)
    x_val = np.arange(index_min, index_max, 1)
    data.index = x_val
    data.to_csv(output_path, header = True, index=True)


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
    return chipnames[0]


def get_labels(start_bin, end_bin, nticks, resolution):
    interval = (end_bin - start_bin)
    ticks = list(np.linspace(0, interval, nticks).astype(int))
    pos = np.linspace(start_bin, end_bin, nticks) * resolution/1000 # kb
    pos = pos.astype(int).astype(str)
    labels = [i + ' kb' for i in pos]
    return ticks, labels


def check_prefix(vector):
    return ["chr" in chrm for chrm in vector]


def check_chrnames(labels_config, labels):
    if np.any(check_prefix(labels)) and not np.all(check_prefix(labels_config)):
        labels_config = ['chr'+ chrm for chrm in labels_config]
    elif not np.any(check_prefix(labels)) and np.all(check_prefix(labels_config)):   
        labels_config = [chrm[3:] for chrm in labels_config]
        
    if set(labels_config).issubset(set(labels)):
        chrnames = labels_config
    elif len(set(labels).intersection(labels_config)) > 0:
        chrnames = np.array(set(labels).intersection(labels_config))
    else:
        self.log.info('Specified choromosomes are not found! Please set correct chromosome names in the configuration file.')
        sys.exit(1)
            
    return chrnames


def check_filenames(hic_files, chipseq_files):
    samplenames = [os.path.split(i)[1].split('.')[0] for i in hic_files]
    hic_files_number = len(hic_files)
    chipseq_files_number = len(chipseq_files)
    if hic_files_number != chipseq_files_number:
        if hic_files_number == 1:
            hic_files = np.repeat(hic_files, chipseq_files_number)
            samplenames = [os.path.split(i)[1].split('.')[0] for i in chipseq_files]
        elif chipseq_files_number == 1:
            chipseq_files = np.repeat(chipseq_files, hic_files_number)
        else:
            print("ERROR: Please provide correct number of input Hi-C and ChIP-seq files")
        
    return np.array(sorted(hic_files)), np.array(sorted(chipseq_files)), np.array(sorted(samplenames))

