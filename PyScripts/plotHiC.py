import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def plotHiC(matrix, segments = None):
    """ Plotting HiC contact matrix
        
        Parameters
        ----------
        ``fname``: array-like type
        HiC matrix
        ``segmentation``: bool
        True or False
        ``segments``: array-like type
        Borders of the TADs
        """
    
    print("Plotting HiC heatmap")
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    sns.heatmap(matrix, cmap = 'Blues')
    
    if segments is not None:
        plt.plot([segments[:, 0], segments[:, 1]], [segments[:, 0], segments[:, 0]], '--', linewidth = 2.5, color = 'firebrick')
        plt.plot([segments[:, 1], segments[:, 1]], [segments[:, 0], segments[:, 1]], '--', linewidth = 2.5, color = 'firebrick')


if __name__ == "__main__":
    correction = sys.argv[1]
    segmets = sys.argv[2]
    
    plotHiC(matrix, segments = None)
