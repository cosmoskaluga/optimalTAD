#!/usr/bin/python
import warnings
warnings.filterwarnings("ignore")

import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import h5py

from scipy import stats
from collections import Counter
from loess.loess_1d import loess_1d

from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression

from pylab import rcParams

from scipy.optimize import curve_fit

plt.rcParams["font.family"] = "serif"


def polynomial_regression(X, y, order):
    poly_reg = PolynomialFeatures(degree = order)
    X_poly = poly_reg.fit_transform(X.reshape(-1, 1))
    
    lin_reg_2 = LinearRegression()
    lin_reg_2.fit(X_poly, y.reshape(-1, 1))
    return lin_reg_2.predict(X_poly)

def fsigmoid(x, a, b, c, d):
    return (1 / (1 + np.exp(a*(x+b))))*c - d


if __name__ == "__main__":
    dirName = sys.argv[1]
    repetition = sys.argv[2]
    chip_seq_filename = "depl_chipseq.merged.binned_short"
    
    ####CHIP-SEQ DATA####
    chrs_drop = ['chr2LHet', 'chr2RHet', 'chr3LHet','chr3RHet','chr4', 'chrU','chrUextra','chrXHet','chrYHet']
    #chrs_drop = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr2LHet', 'chr2RHet', 'chr3LHet', 'chr3RHet', 'chr4', 'chrU', 'chrUextra', 'chrXHet', 'chrYHet']

    #chip_names = ['1','2', '5','9','6','10']
    #index = np.arange(1, 13)
    #index = [str(i) for i in index]
    
    chip_data = pd.read_csv('./ChIP_Seq/' + chip_seq_filename, sep = '\t')
    #chip_data.columns = index
    #chip_data = chip_data[chip_names]

    for chrs in chrs_drop:
        chip_data = chip_data.loc[chip_data['Chr'] != chrs]

    chip_data.index = range(len(chip_data))
    #chip_data.columns = ['Chr', 'Bp', 'Ch1lacz', 'Ch5lacz', 'Ch4tsa', 'Ch8tsa']
    repetition_norm = stats.zscore(chip_data[[repetition]].values)


    if not os.path.exists(dirName):
        os.makedirs('./Output/Score/' + dirName, exist_ok=True)
        os.makedirs('./Output/Pictures/' + dirName, exist_ok=True)

    f = open('./Output/Score/'+dirName+'/output.csv', 'w')
    f.write('Gamma' + '\t' + 'Height_median' + '\t' + 'Height_mean' + '\t' + 'a' + '\t' + 'b' + '\t' + 'c' + '\t' + 'd' + '\n')

    ####TAD DISTANCE######
    for dist_filename in glob.glob('./Output/NumeratedBorders/'+ dirName +'/*.dist'):
        #print('     File: ' + dist_filename[26:])
        dist_data = pd.read_csv(dist_filename, header = None, sep = ',')
        dist_data.columns = ['Chr', 'Bp','Number']

        chrm_index = [1151, 2209, 3437, 4833, 6023 - 68]
        chrm_index = chrm_index[::-1]
    
        for index in chrm_index:
            dist_data = dist_data.drop([index-1])

        for chrs in chrs_drop:
            dist_data = dist_data.loc[dist_data['Chr'] != chrs]

        dist_data.index = range(len(dist_data))


        ####PLOTTING######
        plt.figure(figsize=(10, 7))
        rcParams.update({'font.size': 15})
        plt.style.use('seaborn-darkgrid')

        kb_list = np.arange(-4.0, 4.1, 1.0)
        z = []
        z_interTAD = []
        z_TAD = []
        z_new_mean = []
        for i in kb_list:
            row = dist_data.loc[dist_data['Number'] == i]
            row_indexes = row.index.values.tolist()
            slices = []
            for num in row_indexes:
                if -3 < repetition_norm[num] < 5:
                    #slices.append(repetition_norm[num])
                    slices = slices + list(repetition_norm[num])

            if i <= -1:
                z_interTAD = z_interTAD + slices
            if i >= 1:
                z_TAD = z_TAD + slices

            z.append(np.median(slices))


        #x_pred, y_pred, weigts = loess_1d(kb_list*20, np.array(z), frac=0.55)
        #plt.plot(kb_list*20, y_pred, linewidth = 2.8)
    
        """
        popt, pcov = curve_fit(fsigmoid, kb_list*20, np.array(z), bounds=([0, 0,  0, -1], [3, 10, 3, 1]),  method='dogbox')
        a = popt[0]
        b = popt[1]
        c = popt[2]
        d = popt[3]

        y_points_new = [fsigmoid(i, a, b, c, d) for i in kb_list*20]

        a = str(a)
        b = str(b)
        c = str(c)
        d = str(d)

        y_pred = polynomial_regression(kb_list*20, np.array(z),  3)
        plt.title('Gamma = ' + dist_filename[60:-7], fontsize = 30, y = 1.02)
        plt.plot(kb_list*20, np.array(z), linewidth = 2.8, label = 'reference')
        plt.plot(kb_list*20, y_points_new, linewidth=2.8, label = 'prediction')
        plt.xticks(kb_list*20)
        plt.xlabel('Distance to TAD boundary, kb')
        plt.ylabel('Average Z-score of acetylation values')
        plt.axis([-85, 85, -1.6, 1.6])
        plt.legend()
        plt.savefig('./Output/Pictures/' + dist_filename[26:-5] + '.png', pad_inches = 0.1, dpi = 130)
        """
        
        height_median = np.median(z_interTAD) - np.median(z_TAD)
        height_mean = np.mean(z_interTAD) - np.mean(z_TAD)

        output_gamma = dist_filename[-11:-5]
        if output_gamma[0].isdigit() == False:
            output_gamma = output_gamma[1:]
            if output_gamma[0].isdigit() == False:
                output_gamma = output_gamma[2:]

        ind = 0
        while ind != -1:
            ind = output_gamma.rfind(".", 2, len(output_gamma))
            if ind != -1:
                output_gamma = output_gamma[:ind] + output_gamma[ind+1:]

        a = str(0)
        b = str(0)
        c = str(0)
        d = str(0)
        f.write(output_gamma + '\t' + str(height_median) + '\t' + str(height_mean) + '\t' + a + '\t' + b + '\t' + c + '\t' + d + '\n')
    f.close()





