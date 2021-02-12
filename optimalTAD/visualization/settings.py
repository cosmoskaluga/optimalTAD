import numpy as np
import pandas as pd
from pylab import rcParams
from scipy import ndimage
import matplotlib.pyplot as plt
import logging

from .. optimization import utils

log = logging.getLogger(__name__)

class Plot:
    def __init__(self, path_to_hic, chromosome, start_bin, end_bin, resolution):
        self.path_to_hic = path_to_hic
        self.chromosome = chromosome
        self.start_bin = start_bin
        self.end_bin = end_bin
        self.resolution = resolution

        hic_matrix = np.loadtxt(path_to_hic)
        hic_matrix = hic_matrix[self.start_bin:self.end_bin, self.start_bin:self.end_bin]
        self.matrix = ndimage.rotate(hic_matrix, 45, order = 0, reshape = True, prefilter = False, cval = np.nan)
        
        x_min = self.start_bin
        x_max = self.end_bin
        y_min = self.start_bin
        y_max = np.sqrt(x_max * x_max * 2)
        self.coeff = (y_max - y_min)/(x_max - x_min)
    
    def plotHiC(self, text = 'Hi-C', cmap = 'coolwarm', nticks = 4, figsize = (11, 4)):
        self.position = np.array([0.15, 0.4, 0.8, 0.6])
        p = self.position[1] + (1-self.position[1])/2
        hic_text_pos  = np.array([0.07, p, 0.1, 0.1])
        cbap_pos = np.array([0.9, 0.8, 0.015, 0.22])
        self.fig = plt.figure(figsize=figsize)
        self.nticks = nticks
        
        h_ax = self.fig.add_axes(self.position)
        text_ax = self.fig.add_axes(hic_text_pos)
        cbar_ax = self.fig.add_axes(cbap_pos)
        
        m = self.matrix
        m[m == m.min()] = np.nan
        im = h_ax.imshow(self.matrix, cmap = cmap)
        self.ylim_start = self.matrix.shape[0]//2 - 1
        self.ylim_end = self.ylim_start//2
        h_ax.set(ylim = (self.ylim_start, self.ylim_end))
        h_ax.axis('off')
        
        self.h_ax = h_ax
        
        valmin = np.nanmin(self.matrix)
        valmax = np.nanmax(self.matrix)
        cbar = self.fig.colorbar(im, cax = cbar_ax, ticks = [valmin, valmax], format = '%.3g')
        cbar_ax.tick_params(labelsize = 12)
    
        text_ax.text(0, 0, text, fontsize = 12)
        text_ax.axis('off')

    def plotTAD(self, tad, tad_linewidth = 1.2, tad_linestyle = '--'):
        tad.Start = tad.Start.div(self.resolution)
        tad.End = tad.End.div(self.resolution)
        tad_short = tad.loc[(tad.Start >= self.start_bin) & (tad.End <= self.end_bin)]
        border = tad_short[['Start', 'End']].values
        self.border = border
        
        for sl in border:
            mid_point = (sl[1] - sl[0])/2
            x_val = (np.array([sl[0], sl[0] + mid_point, sl[1]]) - self.start_bin) * self.coeff
            y_val = np.array([self.ylim_start, self.ylim_start - mid_point*self.coeff, self.ylim_start])
            self.h_ax.plot(x_val, y_val, linewidth = tad_linewidth, linestyle = tad_linestyle, color = 'black')


    def plotTrack(self, chip_data,  text = 'ChIP-seq', vline_linewidth = 1., vline_linestyle = 'dashed', fontsize = 12):
        old_pos = self.position[1]
        self.position[1] -= 0.18
        self.position[3] = 0.18
        chip_ax = self.fig.add_axes(self.position)
        p = self.position[1] + (old_pos-self.position[1])/2
        chip_text_pos = np.array([0.07, p, 0.1, 0.1])
        text_ax = self.fig.add_axes(chip_text_pos)

        chip_score = chip_data.Score.values[self.start_bin:self.end_bin]
        x_val = chip_data.Start.values/self.resolution
        x_val = x_val[self.start_bin:self.end_bin]
        chip_ax.plot(x_val, chip_score, color = 'black', zorder = 3)
        chip_ax.fill_between(x_val, chip_score, color = '0.8', zorder = 2)
        chip_ax.set_xlim(self.start_bin - 0.5, self.end_bin)
        chip_ax.axis('off')

        try:
            for sl in self.border.flatten():
                chip_ax.axvline(sl, linewidth = vline_linewidth, linestyle = vline_linestyle, color = 'grey', zorder = 1) 
        except NameError:
            pass

        text_ax.text(0, 0, text, fontsize = fontsize)
        text_ax.axis('off')
    
    def show(self):
        #self.position[1] += 0.04
        self.position[3] = 0.01
        c_ax = self.fig.add_axes(self.position)
        c_ax.tick_params(axis = 'both', bottom = True, left = False,
                         right = False, labelbottom = True, labelleft = False)

        ticks, labels = utils.get_labels(self.start_bin, self.end_bin, self.nticks, self.resolution)
        c_ax.set_xlim(ticks[0], ticks[-1])
        c_ax.set_xticks(ticks)
        c_ax.set_xticklabels(labels, fontsize = 12)
        c_ax.set_ylim(0, 0.02)
        c_ax.set_xlabel(self.chromosome, fontsize = 15)

        self.fig.show()

    def saveplot(self, filename, dpi = 200, bbox_inches = 'tight'):
        self.fig.savefig(filename, dpi = dpi, bbox_inches = bbox_inches)


