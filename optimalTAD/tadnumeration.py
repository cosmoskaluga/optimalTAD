import numpy as np
import pandas as pd
import sys
import os
import glob



class TAD_files:
    def __init__(self, input_path, chrs, gamma_max, stepsize):
        self.input_path = input_path
        self.chrs = chrs
        self.gamma_values = np.linspace(0, gamma_max, round(gamma_max/stepsize) + 1, endpoint=True)
    
    def get_files(self):
        tad_files = {key: [] for key in self.gamma_values}
        for chromosome in self.chrs:
            path = os.path.join(self.input_path, chromosome + '/*')
            filenames = glob.glob(path)
            for gamma, file in zip(self.gamma_values, filenames):
                tad_files[gamma].append(file)
        
        return tad_files



class Distances:
    def __init__(self, resolution):
        self.resolution = resolution
    
    def numeration(self, size, interaction_type = 'tad'):
        if size%2 == 0:
            dist = np.arange(size//2)
            dist = np.concatenate([dist, dist[::-1]])
        else:
            dist = np.arange((size-1)//2)
            dist = np.concatenate([dist, [size//2], dist[::-1]])
        
        if interaction_type == 'intertad':
            dist = (dist + 1) * (-1)
        return dist
    
    def columns(self, first_col, second_col, start, end, interaction_type = 'tad'):
        size = (end-start)/self.resolution
        col1 = np.arange(start, end, self.resolution)
        col2 = self.numeration(int(size), interaction_type)
        first_col = np.append(first_col, col1)
        second_col = np.append(second_col, col2)
        
        return first_col, second_col
    
    def get_distances(self, data, length, chromosome):
        bp = 0
        first_col = np.array([])
        second_col = np.array([])
        for row in data:
            start, end = row[1], row[2] + 1
            if bp != start:
                first_col, second_col = self.columns(first_col, second_col, bp, start, 'intertad')
            first_col, second_col = self.columns(first_col, second_col, start, end, 'tad')
            bp = end
        
        if (length > bp):
            start, end = bp, length
            first_col, second_col = self.columns(first_col, second_col, bp, end, 'intertad')
        
        labels = np.repeat([chromosome], len(first_col))
        out = np.vstack((labels, first_col.astype(int), second_col.astype(int)))
        return np.transpose(out)



class Indexing:
    def __init__(self, chrs, resolution, chr_length, samplename, gamma_max, stepsize):
        self.chr_labels = chrs
        self.resolution = resolution
        self.chr_length = chr_length * self.resolution
        self.input_path = 'output/tads/' + samplename
        self.output_path = 'output/borders/'
        self.distances = Distances(resolution)
        self.__tad_files = TAD_files(self.input_path, self.chr_labels, gamma_max, stepsize).get_files()
    
    def get_indexes(self):
        dict_md = {key: [] for key in self.__tad_files.keys()}
        for gamma in self.__tad_files.keys():
            df = pd.DataFrame()
            lbl_total = self.chr_labels
            for path, length, label in zip(self.__tad_files[gamma], self.chr_length, self.chr_labels):
                data = pd.read_csv(path, header = None, names = ['Chr', 'start', 'end'], sep = '\t')
                if data.empty == True:
                    lbl_total = np.delete(lbl_total, label)
                    print('     There are no TADs for {} chromosome, skipping'.format(label))
                else:
                    data = data[::-1]
                    markdown = self.distances.get_distances(data.values, length, label)
                    df_chr = pd.DataFrame(markdown)
                    df = pd.concat([df, df_chr])
            dict_md[gamma].append(df.values)
            dict_md[gamma].append(lbl_total)
        return dict_md






