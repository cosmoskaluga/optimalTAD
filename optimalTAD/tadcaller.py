import os
import sys
import numpy as np

class armatus:
    def __init__(self, np, resolution, path_to_data, gamma_max, stepsize, output, chromosome):
        self.np = str(np)
        self.resolution = str(resolution)
        self.input = path_to_data
        self.gamma_max = str(gamma_max)
        self.stepsize = str(stepsize)
        self.chromosome = chromosome
        self.output = os.path.join(output, "armatus")
        self.path_to_armatus = './include/ArmatusParallel/src/armatus'
    
    def run(self):
        os.system("mpiexec -np {} {} -r {} -i {} -g {} -o {} -s {} -m -c {}".format(self.np,
                                                                               self.path_to_armatus,
                                                                               self.resolution,
                                                                               self.input,
                                                                               self.gamma_max,
                                                                               self.output,
                                                                               self.stepsize,
                                                                               self.chromosome))

def check_path(path, folder_name, name):
    dirName = os.path.join(path, folder_name, name + "/")
    if not os.path.exists(dirName):
        os.makedirs(dirName, exist_ok=True)
    return dirName


def run_tadcaller(args, data, file):
    path = os.path.join(sys.path[0], "output")
    filename = os.path.split(file)[1]
    
    path_to_hic = check_path(path, "data", filename.split(".")[0])
    
    for chromosome in data.keys():
        path_to_file = os.path.join(path_to_hic, chromosome + "." + args.hic_format)
        submap = data[chromosome]
        np.savetxt(path_to_file, submap, delimiter = '\t', fmt = '%.3f')
        armatus_output = check_path(path, "tads/" + filename.split(".")[0], chromosome)
        caller = armatus(args.np, args.resolution, path_to_file, args.gamma_max, args.stepsize, armatus_output, chromosome).run()






