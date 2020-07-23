import os
import sys
import numpy as np


def armatus(np, resolution, path_to_data, gamma_max, stepsize, output, chromosome):
    path_to_output = os.path.join(output, 'armatus')
    path_to_armatus = './armatus/src/armatus'
    os.system('mpiexec -np {} {} -r {} -i {} -g {} -o {} -s {} -m -c {}'.format(str(np),
                                                                                path_to_armatus,
                                                                                str(resolution),
                                                                                path_to_data,
                                                                                str(gamma_max),
                                                                                path_to_output,
                                                                                str(stepsize),
                                                                                chromosome))








