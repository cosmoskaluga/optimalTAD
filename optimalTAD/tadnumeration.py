import numpy as np
import pandas as pd
import sys
import os
import glob


def get_tad_files(input_path, chrs, gamma_max, stepsize):
    gamma_values = np.linspace(0, gamma_max, round(gamma_max/stepsize) + 1, endpoint=True)
    tad_files = {key: [] for key in gamma_values}
    for chromosome in chrs:
        path = os.path.join(input_path, chromosome + '/*')
        filenames = glob.glob(path)
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
    col2 = numeration(int(size), interaction_type)
    first_col = np.append(first_col, col1)
    second_col = np.append(second_col, col2)
    return first_col, second_col



def get_distances(data, length, chromosome, resolution):
    bp = 0
    first_col = np.array([])
    second_col = np.array([])
    for row in data:
        start, end = row[1], row[2] + 1
        if bp != start:
            first_col, second_col = columns(first_col, second_col, bp, start, resolution, 'intertad')
        first_col, second_col = columns(first_col, second_col, start, end, resolution, 'tad')
        bp = end
        
    if (length > bp):
        start, end = bp, length
        first_col, second_col = columns(first_col, second_col, bp, end, resolution, 'intertad')
        
    labels = np.repeat([chromosome], len(first_col))
    out = np.vstack((labels, first_col.astype(int), second_col.astype(int)))
    return np.transpose(out)



def tadnumeration(chr_labels, resolution, chr_length, samplename, gamma_max, stepsize):
    chr_length = chr_length * resolution
    input_path = 'output/tads/' + samplename
    output_path = 'output/borders/'
    tad_files = get_tad_files(input_path, chr_labels, gamma_max, stepsize)

    dict_md = {key: [] for key in tad_files.keys()}
    for gamma in tad_files.keys():
        df = pd.DataFrame()
        lbl_total = chr_labels
        for path, length, label in zip(tad_files[gamma], chr_length, chr_labels):
            data = pd.read_csv(path, header = None, names = ['Chr', 'start', 'end'], sep = '\t')
            if data.empty == True:
                lbl_total = np.delete(lbl_total, label)
                print('     There are no TADs for {} chromosome, skipping'.format(label))
            else:
                data = data[::-1]
                markdown = get_distances(data.values, length, label, resolution)
                df_chr = pd.DataFrame(markdown)
                df = pd.concat([df, df_chr])
        dict_md[gamma].append(df.values)
        dict_md[gamma].append(lbl_total)
    return dict_md





