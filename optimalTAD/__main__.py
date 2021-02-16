import argparse
import logging
import sys
import time
import glob
import os
from . import logger
from . import config
from . visualization import plot
from . optimization import run

class optimalTAD:
    def __init__(self):
        self.log = logger.initialize_logger()
        self.cfg = config.get_configuration()
        
        parser = argparse.ArgumentParser(description = 'optimalTAD: Topologically Associating Domain optimal set prediction', usage = ''' optimalTAD <command> [<args>]
    
The basic optimalTAD commands are:
    run            Run optimization process
    visualize      Visualize results ''')
        
        parser.add_argument('command', help='Subcommand to run')
        args = parser.parse_args([self.cfg.get('basic', 'mode')])

        arg = sys.argv[1:2]
        if arg:
            if arg[0] in ['run', 'visualize']:
                args.command = arg[0]
        
        if args.command not in ['run', 'visualize']:
            self.log.info('Unrecognized command!')
            parser.print_help()
            sys.exit(1)

        getattr(self, args.command)()

    def run(self):
        start_time = time.time()

        hicpath = self.cfg.get('basic','hic')
        chippath = self.cfg.get('basic','chipseq')
        hicpath = glob.glob(os.path.expanduser(hicpath))
        chippath = glob.glob(os.path.expanduser(chippath))
        
        parser = argparse.ArgumentParser(description='Run optimization process')
        parser.add_argument('--hic', type = str, nargs='+', default = sorted(hicpath), help = 'Iteratively corrected Hi-C data')
        parser.add_argument('--chipseq', type = str, nargs = '+', default = sorted(chippath), help = 'ChIP-seq data')
        parser.add_argument('--np', type = int, default = int(self.cfg['basic']['np']), help = 'Number of processors')
        parser.add_argument('--resolution', type = int, default = int(self.cfg['basic']['resolution']), help = 'Resolution')
        parser.add_argument('--stepsize', type = float, default = float(self.cfg['basic']['stepsize']), help = 'Step size to increment gamma parameter')
        parser.add_argument('--gamma_max', type = float, default = float(self.cfg['basic']['gamma_max']), help = 'Max gamma parameter')
        parser.add_argument('--hic_format', type = str, default = self.cfg['basic']['hic_format'], help = 'Hi-C matrices input format for armatus')
        parser.add_argument('--empty_row_imputation', action = 'store_true', help = 'Missing rows (and columns) imputation')
        parser.add_argument('--truncation',  action = 'store_true', help = 'Value truncation of Hi-C-matrix')
        parser.add_argument('--log2_hic', action = 'store_true', help = 'log2 transformation of Hi-C matrix')
        parser.add_argument('--log2_chip', action = 'store_true', help = 'log2 transformation of ChIP-Seq track')
        parser.add_argument('--zscore_chip', action = 'store_true', help = 'Z-score transformation of ChIP-Seq track')
        parser.set_defaults(empty_row_imputation = eval(self.cfg['basic']['empty_row_imputation']))
        parser.set_defaults(truncation = eval(self.cfg['basic']['truncation']))
        parser.set_defaults(log2_hic = eval(self.cfg['basic']['log2_hic']))
        parser.set_defaults(log2_chip = eval(self.cfg['basic']['log2_chip']))
        parser.set_defaults(zscore_chip = eval(self.cfg['basic']['zscore_chip']))
        args = parser.parse_args(sys.argv[2:])
        run.main(args, self.cfg, self.log)
        
        cpu_time = round(time.time()-start_time, 2)
        self.log.info('Execution time: {} sec'.format(cpu_time))
    
    def visualize(self):
        start_time = time.time()
        parser = argparse.ArgumentParser(description='Visualize results')
        parser.add_argument('--samplename', type = str, help = 'Samplename of Hi-C data')
        parser.add_argument('--region', type = str, help = 'Genome region')
        parser.add_argument('--resolution', type = int, help = 'Resolution')
        parser.add_argument('--chipseq', type = str, default = False, help = 'ChIP-seq data')
        parser.add_argument('--log2_chip', action = 'store_true', help = 'log2 transformation of ChIP-Seq track')
        parser.add_argument('--zscore_chip', action = 'store_true', help = 'Z-score transformation of ChIP-Seq track')
        parser.add_argument('--rnaseq', type = str, default = False, help = 'RNA-seq data')
        parser.set_defaults(log2_chip = False)
        parser.set_defaults(zscore_chip = False)
        args = parser.parse_args(sys.argv[2:])
        #plot.main(args, self.cfg['visualization'], self.log)
        
        cpu_time = round(time.time()-start_time, 2)
        self.log.info('Execution time: {} sec'.format(cpu_time))


    
if __name__ == '__main__':
    optimalTAD()
