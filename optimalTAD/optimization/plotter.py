import os
import numpy as np
import pandas as pd
from pylab import rcParams
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import ScalarFormatter, AutoMinorLocator


def stylize_axes(ax):
    """ 
        Mtplotlib axes settings  

        Parameters 
        ----------
        ``ax`` : matplotlib-like object
            An object containing the matplotlib figure 
    """
    ax.yaxis.major.formatter._useMathText = True
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis='both', which='major', labelsize=15)


def plotAmplitude(data, output_path, dpi = 200):
    """ 
        Plotting a difference in median ChIP-seq value between inter-TADs and TADs for each value of the optimized parameter 

        Parameters 
        ----------
        ``data`` : dataframe
            A dataframe of delta h differences in median ChIP-seq signal per gamma (armatus) or window size (IS)
        ``output_path`` : int
            A path to save figure 
        ``dpi`` : int
            Figure resolution 
        
        Returns
        -------
        ``ax`` : matplotlib-like object
            An object containing the generated amplitude figure 
    """
    sns.set_palette(sns.color_palette('Set1'))
    ci_df = data.filter(like='CI', axis=1)
    sample_df = data.drop(list(data.filter(like='CI', axis=1).columns), axis=1)
    x_val = data.Gamma
    
    fig, ax = plt.subplots(figsize=(9, 7))
    color_list = list(sns.color_palette("pastel", n_colors = len(sample_df.columns[1:])).as_hex())
    for sample, color in zip(sample_df.columns[1:], color_list):
        y_val = data[sample].values
        ax.plot(x_val, y_val, linewidth = 2.5, ls = '--', color = color, marker = 'o', label = sample)
        ci_df_sample = ci_df.filter(like=sample, axis=1)
        ax.fill_between(x_val, ci_df_sample.iloc[:,0], ci_df_sample.iloc[:,1], color = color, alpha = .15)

    ax.set_xlabel('Gamma', fontsize = 20)
    ax.set_ylabel('Amplitude', fontsize = 20)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, fontsize = 12)

    stylize_axes(ax)
    ax.grid(linestyle=':', linewidth='0.2', color='black')
    
    if output_path:
        path = os.path.dirname(output_path)
        if not os.path.exists(path) and path != '':
            os.makedirs(path, exist_ok=True)
        ax.figure.savefig(output_path, dpi=dpi, bbox_inches='tight')

    return ax


def plotStair(stair_df, best_gamma, output_path, index_min = -5, index_max = 5, dpi = 200, path_to_stair_dataframe = None):
    """ 
        Plotting median ChIP-seq values per each distance to a TAD boundary

        Parameters 
        ----------
        ``stair_df`` : dataframe
            Median ChIP-seq values per each index used (per each distance to the nearest TAD boundary)
        ``best_gamma`` : dataframe
            Optimized value of the gamma (armatus) or window size (IS) parameter
        ``index_min`` : int
            Minimal distance value to consider while calculating median of ChIP-seq signal
        ``index_max`` : ind
            Maximal distance value to consider while calculating median of ChIP-seq signal
        ``output_path`` : str
            A path to save figure 
        ``dpi`` : int
            Figure resolution 
        ``path_to_stair_dataframe`` : str
            A path to save the dataframe 
        
        Returns
        -------
        ``ax`` : matplotlib-like object
            An object containing the generated stair figure 
    """
    sns.set_palette(sns.color_palette('Set1'))
    fig, ax = plt.subplots(figsize=(9, 7))
    color_list = list(sns.color_palette("pastel", n_colors = len(list(stair_df.keys()))).as_hex())
    for name, gamma, color in zip(list(stair_df.keys()), best_gamma, color_list):
        x_val = np.arange(index_min, index_max, 1)
        y_val = stair_df[name].values
    
        if np.isnan(y_val).any() == True:
            idx = np.argwhere(np.isnan(y_val))
            idx = np.concatenate(idx)
            x_val = np.delete(x_val, idx)
            y_val = np.delete(y_val, idx)
    
        poly = np.polyfit(x_val, y_val, 3)
        poly_y = np.poly1d(poly)(x_val)
        ax.plot(x_val, poly_y, linewidth = 2.5, color = color, label = name + ' ($\gamma$ = ' + str(gamma) + ')')

    ax.set_xlabel('Distance to TAD boundary, kb', fontsize = 20)
    ax.set_ylabel('Median z-score of acetylation values', fontsize = 20)
    ax.grid(linestyle=':', linewidth='0.2', color='black')
    
    stylize_axes(ax)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, fontsize = 12)

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




