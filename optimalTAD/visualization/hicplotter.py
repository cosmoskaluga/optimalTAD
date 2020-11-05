import numpy as np
import pandas as pd
from pylab import rcParams
from scipy import ndimage
import matplotlib.pyplot as plt
import logging


from .. optimization import chipseqloader


class Plot:
    def __init__(self, path_to_hic, region, resolution, log2_chip, zscore_chip, path_to_chipseq = False):
        self.path_to_hic = path_to_hic
        self.resolution = resolution
        self.path_to_chipseq = path_to_chipseq
        self.log2_chip = log2_chip
        self.zscore_chip = zscore_chip
        
        hic_matrix = np.loadtxt(path_to_hic)
        
        region_split = region.split(':')
        try:
            self.chromosome = region_split[0]
            if len(region_split) == 2:
                coordinates = region_split[1].split('-')
                self.start_bin = int(coordinates[0].replace(',',''))
                self.end_bin = int(coordinates[1].replace(',',''))
            
                self.start_bin = int(self.start_bin/self.resolution)
                self.end_bin = int(self.end_bin/self.resolution)
            else:
                self.start_bin = 0
                self.end_bin = np.shape(hic_matrix)[0]
        except ValueError:
            print('Invalid format of chromosome coordinates!')
        
        hic_matrix = hic_matrix[self.start_bin:self.end_bin, self.start_bin:self.end_bin]
        self.matrix = ndimage.rotate(hic_matrix, 45, order=0, reshape=True, prefilter=False, cval=np.nan)
        
        if not path_to_chipseq == False:
            ChipSeqLoader = chipseqloader.ChipSeq(path_to_chipseq)
            chip_data = ChipSeqLoader(self.log2_chip, self.chromosome, self.zscore_chip)
            self.chip_data = chip_data.loc[chip_data.Chr.isin([self.chromosome])]
        
        x_min = self.start_bin
        x_max = self.end_bin
        y_min = self.start_bin
        y_max = np.sqrt(x_max*x_max*2)
        self.coeff = (y_max - y_min)/(x_max - x_min)
    
    def plotHiC(self, figsize = (11, 4), text = 'Hi-C', cmap = 'coolwarm', nticks = 4):
        #coolwarm
        cmap = 'coolwarm'
        self.fig = plt.figure(figsize=figsize)
        heatmap_pos=[0.15, 0.4, 0.8, 0.7]
        chrom_pos=[0.15, 0.14, 0.8, 0.010]
        text_pos = [0.07, 0.75, 0.1, 0.1]
        cbap_pos = [0.9, 0.8, 0.015, 0.22]
        
        h_ax = self.fig.add_axes(heatmap_pos)
        c_ax = self.fig.add_axes(chrom_pos)
        text_ax = self.fig.add_axes(text_pos)
        cbar_ax = self.fig.add_axes(cbap_pos)
        
        m = self.matrix
        m[m == m.min()] = np.nan
        im = h_ax.imshow(self.matrix, cmap = cmap)
        self.ylim_start = self.matrix.shape[0]//2 - 1
        self.ylim_end = self.ylim_start//2
        h_ax.set(ylim = (self.ylim_start, self.ylim_end))
        h_ax.axis('off')
        
        c_ax.tick_params(axis='both', bottom=True, top=False, left=False,
                         right=False, labelbottom=True, labeltop=False,
                         labelleft=False, labelright=False)
            
        interval = (self.end_bin - self.start_bin)
        ticks = list(np.linspace(0, interval, nticks).astype(int))
        pos = np.linspace(self.start_bin, self.end_bin, nticks) * self.resolution/1000
        pos = pos.astype(int).astype(str)
        labels = [i + ' kb' for i in pos]
        
        c_ax.set_xlim(ticks[0], ticks[-1])
        c_ax.set_xticks(ticks)
        c_ax.set_xticklabels(labels, fontsize=12)
        
        c_ax.set_ylim(0, 0.02)
        c_ax.set_xlabel(self.chromosome, fontsize = 15)
        self.h_ax = h_ax
        
        valmin = np.nanmin(self.matrix)
        valmax = np.nanmax(self.matrix)
        cbar = self.fig.colorbar(im, cax=cbar_ax, ticks=[valmin, valmax], format='%.3g')
        cbar_ax.tick_params(labelsize=12)
    
        text_ax.text(0, 0, text, fontsize = 12)
        text_ax.axis('off')

    def plotTAD(self, path_to_tad, vline_linewidth = 1., vline_linestyle = 'dashed', tad_linewidth = 1.2, tad_linestyle = '--'):
        tad = pd.read_csv(path_to_tad, header = None, names = ['Chr', 'Start', 'End'], sep = '\t')
        tad = tad[::-1]
        tad.End = tad.End + 1
        tad.Start = tad.Start.div(self.resolution)
        tad.End = tad.End.div(self.resolution)
        self.tad = tad
        
        tad_short = self.tad.loc[(self.tad.Start >= self.start_bin) & (self.tad.End <= self.end_bin)]
        border = tad_short[['Start', 'End']].values
        
        for sl in border:
            mid_point = (sl[1] - sl[0])/2
            x_val = (np.array([sl[0], sl[0] + mid_point, sl[1]]) - self.start_bin) * self.coeff
            y_val = np.array([self.ylim_start, self.ylim_start - mid_point*self.coeff, self.ylim_start])
            self.h_ax.plot(x_val, y_val, linewidth = tad_linewidth, linestyle = tad_linestyle, color = 'black')
            self.chip_ax.axvline(sl[0], linewidth = vline_linewidth, linestyle = vline_linestyle, color = 'grey', zorder = 1)
            self.chip_ax.axvline(sl[1], linewidth = vline_linewidth, linestyle = vline_linestyle, color = 'grey', zorder = 1)

    def plotChiPSeqTrack(self, text = 'ChIP-seq', fontsize = 12):
        chip_pos = [0.15, 0.15, 0.8, 0.32]
        text_pos = [0.05, 0.25, 0.1, 0.1]
        chip_ax = self.fig.add_axes(chip_pos)
        text_ax = self.fig.add_axes(text_pos)
        
        chip_score = self.chip_data.Score.values[self.start_bin:self.end_bin]
        x_val = self.chip_data.Start.values/self.resolution
        x_val = x_val[self.start_bin:self.end_bin]
        chip_ax.plot(x_val, chip_score, color = 'black', zorder = 3)
        chip_ax.fill_between(x_val, chip_score, color = '0.8', zorder = 2)
        chip_ax.set_xlim(self.start_bin - 0.5, self.end_bin)
        self.chip_ax = chip_ax
        chip_ax.axis('off')
    
        text_ax.text(0, 0, text, fontsize = fontsize)
        text_ax.axis('off')
    
    def saveplot(self, filename, dpi = 200, bbox_inches = 'tight'):
        self.fig.savefig(filename, dpi = dpi, bbox_inches = bbox_inches)
    
    def show(self):
        self.fig.show()


