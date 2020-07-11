import pandas as pd
import os
import glob
import logging

log = logging.getLogger(__name__)


def optimal_gamma(data):
    samples = data.columns[1:]
    gamma = data.Gamma
    for sample in samples:
        maxindex = data[sample].idxmax()
    return round(gamma[maxindex], 2)


def progressbar (iteration, total):
    prefix = " "*21
    suffix = "complete"
    length = 47
    decimals = 1
    fill = '█'
    
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    if iteration == total:
        print()

def check_path(path, folder_name, name = None):
    dirName = os.path.join(path, folder_name, name + "/")
    if not os.path.exists(dirName):
        os.makedirs(dirName, exist_ok=True)
    else:
        for old_file in glob.glob(os.path.join(dirName, "*")):
            os.remove(old_file)
    return dirName
