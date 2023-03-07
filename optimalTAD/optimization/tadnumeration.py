import numpy as np
import pandas as pd
import sys
import os
import glob
import logging

log = logging.getLogger(__name__)


def sort_files(arr):
    arr = np.array(arr)
    gamma_list = np.array([])
    for i in arr:
        gamma_str = i.split('gamma.')[1]
        gamma_str = gamma_str.split('.txt')[0]
        gamma_list = np.append(gamma_list, gamma_str)
    return arr[np.argsort(gamma_list)]



def get_tad_files(input_path, chrs, gamma_max, stepsize):
    gamma_values = np.linspace(0, gamma_max, round(gamma_max/stepsize) + 1, endpoint = True).round(3)
    tad_files = {key: [] for key in gamma_values}
    for chromosome in chrs:
        path = os.path.join(input_path, chromosome + '/*')
        #filenames = glob.glob(path)
        filenames = sort_files(glob.glob(path))
        for gamma, file in zip(gamma_values, filenames):
            tad_files[gamma].append(file)
        
    return tad_files


    
def numeration(size, interaction_type = 'tad'):
    if size%2 == 0:
        dist = np.arange(size//2)
        dist = np.concatenate([dist, dist[::-1]])
    else:
        dist = np.arange((size-1)//2)
        dist = np.concatenate([dist, [size//2], dist[::-1]])
        
    if interaction_type == 'intertad':
        dist = (dist + 1) * (-1)
    return dist



def columns(first_col, second_col, start, end, resolution, interaction_type = 'tad'):
    size = (end-start)/resolution
    col1 = np.arange(start, end, resolution)
    
    if interaction_type != 'skip':
        col2 = numeration(int(size), interaction_type).astype(str)
    else:
        col2 = np.repeat(np.nan, size)

    first_col = np.append(first_col, col1)
    second_col = np.append(second_col, col2)
    return first_col, second_col.astype(str)



def get_distances(data, length, chromosome, resolution):
    first_col = np.array([])
    second_col = np.array([])
    
    bp = 0
    start_bin = int(data[0][1])
    if start_bin != 0:
        first_col, second_col = columns(first_col, second_col, bp, start_bin, resolution, interaction_type = 'skip')
    bp = start_bin

    for row in data:
        start_bin, end_bin = row[1], row[2] + 1
        if bp != start_bin:
            first_col, second_col = columns(first_col, second_col, bp, start_bin, resolution, interaction_type = 'intertad')
        first_col, second_col = columns(first_col, second_col, start_bin, end_bin, resolution, interaction_type = 'tad')
        bp = end_bin

    if (length > bp):
        start, end = bp, length
        first_col, second_col = columns(first_col, second_col, bp, end, resolution, interaction_type = 'skip')
    
    labels = np.repeat([chromosome], len(first_col))
    out = np.vstack((labels, first_col.astype(int), second_col))
    return np.transpose(out)



def get_numeration(chr_labels, resolution, chr_length, samplename, gamma_max, stepsize):
    chr_length = chr_length * resolution
    input_path = os.path.join(os.path.realpath('.'), 'output/tads/') + samplename
    output_path = os.path.join(os.path.realpath('.'),'output/borders/')
    tad_files = get_tad_files(input_path, chr_labels, gamma_max, stepsize)

    dict_md = {key: [] for key in tad_files.keys()}
    chrs_counter = 0
    for gamma in tad_files.keys():
        df = pd.DataFrame()
        lbl_total = chr_labels
        for path, length, label in zip(tad_files[gamma], chr_length, chr_labels):
            data = pd.read_csv(path, header = None, names = ['Chr', 'Start', 'End'], sep = '\t')
            if data.empty == True:
                index = np.argwhere(lbl_total == label)
                lbl_total = np.delete(lbl_total, index)
                #print('     There are no TADs for {} chromosome, skipping'.format(label))
            else:
                data = data[::-1]
                markdown = get_distances(data.values, length, label, resolution)
                df_chr = pd.DataFrame(markdown)
                df = pd.concat([df, df_chr])
                chrs_counter += 1
        dict_md[gamma].append(df.values)
        dict_md[gamma].append(lbl_total)

    if chrs_counter == 0:
        log.error('There are no TADs identified across resolutions and replicates of Hi-C data!')
        sys.exit(1)

    return dict_md, tad_files





