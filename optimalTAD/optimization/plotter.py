import os
import numpy as np
import pandas as pd
from pylab import rcParams
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import ScalarFormatter, AutoMinorLocator


def stylize_axes(ax):
    ax.yaxis.major.formatter._useMathText = True
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis='both', which='major', labelsize=15)


def plotAmplitude(data, output_path = 'amplitude.png', dpi = 200):
    sns.set_palette(sns.color_palette('Set1'))
    samples = data.columns[1:]
    x_val = data.Gamma
    
    fig, ax = plt.subplots(figsize=(9, 7))
    for sample in samples:
        y_val = data[sample].values
        ax.plot(x_val, y_val, linewidth = 2.5, ls='--', marker='o', label = sample)

    ax.set_xlabel('Gamma', fontsize = 20)
    ax.set_ylabel('Amplitude', fontsize = 20)
    ax.legend(frameon=False, fontsize = 12)
    stylize_axes(ax)
    ax.grid(linestyle=':', linewidth='0.3', color='black')
    
    if output_path:
        path = os.path.dirname(output_path)
        if not os.path.exists(path) and path != '':
            os.makedirs(path, exist_ok=True)
        ax.figure.savefig(output_path, dpi=dpi, bbox_inches='tight')

    return ax


def plotStair(stair_df, best_gamma, index_min = -5, index_max = 5, output_path = 'stair.png', dpi = 200, path_to_stair_dataframe = None):
    sns.set_palette(sns.color_palette('Set1'))
    fig, ax = plt.subplots(figsize=(9, 7))
    for name, gamma in zip(list(stair_df.keys()), best_gamma):
        x_val = np.arange(index_min, index_max, 1)
        y_val = stair_df[name].values
    
        if np.isnan(y_val).any() == True:
            idx = np.argwhere(np.isnan(y_val))
            idx = np.concatenate(idx)
            x_val = np.delete(x_val, idx)
            y_val = np.delete(y_val, idx)
    
        poly = np.polyfit(x_val, y_val, 3)
        poly_y = np.poly1d(poly)(x_val)
        ax.plot(x_val, poly_y, linewidth = 2.5, label = name + ' ($\gamma$ = ' + str(gamma) + ')')

    ax.set_xlabel('Distance to TAD boundary, kb', fontsize = 20)
    ax.set_ylabel('Median z-score of acetylation values', fontsize = 20)
    ax.grid(linestyle=':', linewidth='0.3', color='black')
    
    stylize_axes(ax)
    ax.legend(frameon = False, fontsize = 12)

    if output_path:
        path = os.path.dirname(output_path)
        if not os.path.exists(path) and path != '':
            os.makedirs(path, exist_ok=True)
        ax.figure.savefig(output_path, dpi=dpi, bbox_inches='tight')

    if path_to_stair_dataframe:
        x_val = np.arange(index_min, index_max, 1)
        stair_df.index = x_val
        stair_df.to_csv(path_to_stair_dataframe, header = True, index=True)

    return ax




