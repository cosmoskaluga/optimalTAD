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


def plotAmplitude(data, output_path = None, dpi = 200):
    sns.set_palette(sns.color_palette("Set1"))
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
    
    #plt.show()
    
    if output_path:
        path = os.path.dirname(output_path)
        if not os.path.exists(path) and path != '':
            os.makedirs(path, exist_ok=True)
        ax.figure.savefig(output_path, dpi=dpi, bbox_inches='tight')

    return ax




