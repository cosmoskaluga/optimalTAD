import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import glob
from pylab import rcParams

def gamma_vs_amplitude(folder_names):
    plt.figure(figsize=(10, 7))
    rcParams.update({'font.size': 15})
    plt.style.use('seaborn-darkgrid')
    
    for folder in folder_names:
        out = pd.read_csv(folder + '/output.csv', sep = '\t')
        h = out['Height_median'].values
        gamma = out['Gamma'].values
        plt.plot(gamma, h, linewidth = 2.5,  label = folder[15:])

    plt.xlabel('Gamma', labelpad = 6)
    plt.ylabel('Amplitude of Z-score of average \n acetylation values', fontsize = 17, labelpad = 4)
    plt.legend()
    plt.savefig('./Output/Pictures/gamma_vs_amplitude.png', dpi = 400)

def sigmoidal_order(folder_names):
    plt.figure(figsize=(14, 7))
    rcParams.update({'font.size': 15})
    plt.style.use('seaborn-darkgrid')
    for folder in folder_names:
        out = pd.read_csv(folder + '/output.csv', sep = '\t')
        plt.plot(out['Gamma'].values, out[['a']].values, '-o', linewidth = 2.5)

    plt.xticks(range(0, 6))
    plt.xlabel('Gamma', fontsize = 25, labelpad = 6)
    plt.ylabel('$Alpha$', fontsize = 25, labelpad = 4)
    #plt.legend()
    plt.savefig('./Output/Pictures/sigmoidal_order.png', dpi = 400)

def best_gamma(folder_names):
    for folder in folder_names:
        out = pd.read_csv(folder + '/output.csv', sep = '\t')
        out = out.sort_values(by=['Height_median'], ascending=False)
    
        slope_cutoff = 0.05

        #slope_factors = out[['a']].values
        gamma_values = out[['Gamma']].values
        #factor = slope_factors[0][0]
        c = 0

        #while factor < slope_cutoff:
        #    if c == len(slope_factors):
        #        break
        #    else:
        #        c+=1
        #        factor = slope_factors[c][0]

        gamma_best = gamma_values[c][0]
        print('     Best gamma for ' + folder[15:] + ' is '+ str(gamma_best))

if __name__ == "__main__":
    plt.rcParams["font.family"] = "serif"

    folder_names = glob.glob('./Output/Score/*')
    gamma_vs_amplitude(folder_names)
    sigmoidal_order(folder_names)
    best_gamma(folder_names)
    
