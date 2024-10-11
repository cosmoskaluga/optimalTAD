import pandas as pd
import numpy as np
import os
import glob
import shutil
import logging
import sys

log = logging.getLogger(__name__)


def optimal_gamma(data):
    samples = data.columns[1]
    gamma = data.Gamma
    maxindex = data[samples].idxmax()
    
    return round(gamma[maxindex], 2)


def select_optimal_tads(tads, optimal_gamma, samplename, output_path, mammal = False):
    dst = os.path.join(check_path(output_path, '', 'optimal_gamma'), samplename + '.domains.tad')
    if not mammal:
        merged_tads = pd.DataFrame(columns = ['Chr', 'Start', 'End'])
        for src in tads[optimal_gamma]:
            tad = pd.read_csv(src, header = None, names = ['Chr', 'Start', 'End'], sep = '\t')
            tad = tad[::-1]
            tad.End += 1
            merged_tads = pd.concat([merged_tads, tad])
    else:
        merged_tads = tads[optimal_gamma]

    merged_tads.to_csv(dst, header = True, index = False)


def save_stairs(data, index_min, index_max, output_path):
    data = pd.DataFrame(data)
    x_val = np.arange(index_min, index_max, 1)
    data.index = x_val
    data.to_csv(output_path, header = True, index=True, float_format='%.5f')


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

def uniquify_path(path):
    """
    Make unique filepath by adding an integer suffix.
    Borrowed from: https://stackoverflow.com/questions/13852700/create-file-but-if-name-exists-add-number
    """
    filename, extension = os.path.splitext(path)
    counter = 1

    while os.path.exists(path):
        path = filename + "-" + str(counter) + extension
        counter += 1

    return path

def check_path(path, folder_name, name = None):
    dirName = os.path.join(path, folder_name, name + "/")

    if not os.path.exists(dirName):
        os.makedirs(dirName, exist_ok=True)

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
        
    if len(set(labels) & set(labels_config)) > 0:
        chrnames = list(set(labels) & set(labels_config))
    else:
        log.info('ERROR: Specified choromosomes are not found! Please set correct chromosome names in the configuration file.')
        sys.exit(1)

    return chrnames


def check_filenames(hic_files, chipseq_files):
    samplenames = [os.path.splitext(os.path.split(i)[1])[0] for i in hic_files]
    hic_files_number = len(hic_files)
    chipseq_files_number = len(chipseq_files)
    if hic_files_number != chipseq_files_number:
        if hic_files_number == 1:
            hic_files = np.repeat(hic_files, chipseq_files_number)
            samplenames = [os.path.splitext(os.path.split(i)[1])[0] for i in chipseq_files]
        elif chipseq_files_number == 1:
            chipseq_files = np.repeat(chipseq_files, hic_files_number)
        else:
            log.info("ERROR: Please provide a correct number of input Hi-C and ChIP-seq files")
            sys.exit(1)
        
    return np.array(sorted(hic_files)), np.array(sorted(chipseq_files)), np.array(sorted(samplenames))


def save_table(table, output, output_category, name):
    path_to_table = os.path.join(check_path(output, '', output_category), name + ".csv")
    table.to_csv(path_to_table, header = True, index=False, float_format='%.5f')