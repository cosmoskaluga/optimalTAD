import os
import numpy as np
import pandas as pd
from pylab import rcParams
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter,AutoMinorLocator


def stylize_axes(ax):
    ax.yaxis.major.formatter._useMathText = True
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis='both', which='major', labelsize=15)


def plotAmplitude(data, samplename = None, output_path = None, dpi = None):
    x_val = list(data.keys())
    y_val = list(data.values())
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.plot(x_val, y_val, linewidth = 2.5, ls='--', marker='o',  label = samplename)
    ax.set_xlabel('Gamma', fontsize = 20)
    ax.set_ylabel('Amplitude', fontsize = 20)
    stylize_axes(ax)
    ax.legend(frameon=False, loc='upper right', fontsize = 12)
    
    plt.show()
    
    if output_path:
        path = os.path.dirname(output_path)
        if not os.path.exists(path):
            os.makedirs(path, exist_ok=True)
        ax.figure.savefig(output_path, dpi=dpi, bbox_inches='tight')

    return ax




