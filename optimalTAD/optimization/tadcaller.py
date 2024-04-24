import os
import sys
import numpy as np
import bioframe
from . import utils



def armatus(np, resolution, path_to_data, gamma_max, stepsize, output, chromosome):
    path_to_output = os.path.join(output, 'armatus')
    #path_to_armatus = 'armatus'
    path_to_armatus = os.path.join(os.path.dirname( __file__ ), '..', 'armatus')
    os.system('mpiexec -n {} {} -r {} -i {} -g {} -o {} -s {} -m -c {}'.format(str(np),
                                                                                path_to_armatus,
                                                                                str(resolution),
                                                                                path_to_data,
                                                                                str(gamma_max),
                                                                                path_to_output,
                                                                                str(stepsize),
                                                                                chromosome))



def run_armatus(args, chromsize, samplename):
    path = os.path.join(os.path.realpath('.'), 'output')
    path_to_hic = os.path.join(path, 'data', samplename + '/')
    num = 0
    for chromosome in chromsize.keys():
        utils.progressbar(num, len(chromsize.keys()))
        path_to_file = os.path.join(path_to_hic, chromosome + '.' + args.hic_format)
        armatus_output = utils.check_path(path, 'tads/' + samplename, chromosome)
        caller = armatus(args.np,
                                   args.resolution,
                                   path_to_file,
                                   args.gamma_max,
                                   args.stepsize,
                                   armatus_output,
                                   chromosome)
        num += 1
    utils.progressbar(num, len(chromsize.keys()))



def ins_table2tads(ins):
    list_of_sizes = [i for i in ins.columns if i.find('is_boundary') == 0]

    for size in list_of_sizes:
        ins.loc[(ins[size] == True) & (ins[f"is_bad_bin"] == True), size] = False

    blacklist = ins.loc[ins.is_bad_bin == True, ["chrom", "start", "end"]]
    blacklist.columns = ['Chr', 'Start', 'End']
    window_sizes = [i.split("_")[2] for i in ins.columns if i.find('is_boundary') == 0]

    tad_files = dict()

    for window_size in window_sizes:
        ins_ws = ins[ins[f"is_boundary_{window_size}"] == False]
        tads = bioframe.merge(ins_ws)
        tads = tads[(tads["end"] - tads["start"]) <= 2000000].reset_index(drop=True)
        tads = tads[['chrom', 'start', 'end']]
        tads.columns = ['Chr', 'Start', 'End']
        tad_files[window_size] = tads

    return tad_files, blacklist 



def run_IS(path, args, set_chromosomes, ignore_diags = None, clr_weight_name = "weight", min_frac_valid_pixels = 0.66, min_dist_bad_bin = 0, verbose = True):
    from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
    import warnings
    warnings.simplefilter('ignore', category = NumbaDeprecationWarning)
    warnings.simplefilter('ignore', category = NumbaPendingDeprecationWarning)
    import cooltools
    import cooler

    output_path = os.path.join(os.path.realpath('.'), 'output')
    utils.check_path('', '', output_path)

    extension = os.path.splitext(path)[1]
    if extension == '.mcool':
        suffix = '::resolutions/' + str(args.resolution)
        path += suffix

    clr = cooler.Cooler(path) 
    windows = (np.arange(args.window_size_min, args.window_size_max, args.resolution)).astype(int)
    insulation_table = cooltools.insulation(clr,    
                                            window_bp = windows, 
                                            nproc = args.np, 
                                            ignore_diags = ignore_diags, 
                                            clr_weight_name = clr_weight_name, 
                                            min_frac_valid_pixels = min_frac_valid_pixels, 
                                            verbose = verbose)

    if set_chromosomes == 'None':
        labels = clr.chromnames
    else:
        labels_config = set_chromosomes.split(',')
        labels = utils.check_chrnames(labels_config, clr.chromnames)
    
    insulation_table = insulation_table.loc[insulation_table['chrom'].isin(labels)]

    return ins_table2tads(insulation_table)








