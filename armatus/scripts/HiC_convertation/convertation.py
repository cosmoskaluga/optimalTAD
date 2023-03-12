#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import itertools
import gzip
import h5py
import warnings

warnings.filterwarnings("ignore")

def plotting_data(arr, name):
    print("Plotting HiC heatmap for " + name)
    fig = plt.figure(figsize=(14, 14))
    ax = fig.add_subplot(111)
    im = ax.matshow(arr, cmap='YlOrRd')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.2)
    fig.colorbar(im, cax=cax)
    plt.savefig('HiC_map_' + name + '.png', pad = 0.1, dpi = 200)

def output_data(name, arr, output_format):
    print("Saving"+" ." + output_format + " file for " + name)
    if output_format == "txt":
        np.savetxt('chromosome/chr' + name + '.txt', arr, fmt = '%0.3e', delimiter='\t')
    else:
        np.savetxt('chromosome/chr' + name + '.txt.gz', arr, fmt = '%0.3e', delimiter='\t')


filename = sys.argv[1]
output_format = sys.argv[2]

#file_names = ['IC-heatmap-20K_(S2-20).hdf5']

chromosome = ['0 0', '1 1', '2 2', '3 3', '4 4', '5 5']
chromosome_names = ['2L', '2R', '3L', '3R', '4', 'X']

f = h5py.File(filename, 'r')
for chrm, name in zip(chromosome, chromosome_names):
    arr = f[chrm][()]
    arr[arr < 0.5] = 0.5
    arr[arr > 1024] = 1024
    arr = np.log2(arr)
    output_data(name, arr, output_format)
    plotting_data(arr, name)


