import argparse
import logging
import sys
import time
from . import calculate
from . import plotting
from . import logger
from . import config

class optimalTAD:
    def __init__(self):
        
        self.log = logger.initialize_logger()
        self.cfg = config.get_configuration()
        
        parser = argparse.ArgumentParser(description = 'optimalTAD: Topologically Associating Domain optimal set prediction', usage=''' optimalTAD <command> [<args>]
    
The basic optimalTAD commands are:
    run            Run optimization process
    visualize      Visualize results ''')
        
        parser.add_argument('command', help='Subcommand to run')
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            self.log.info('Unrecognized command!')
            parser.print_help()
            sys.exit(1)
        getattr(self, args.command)()

    def run(self):
        start_time = time.time()
        parser = argparse.ArgumentParser(description='Run optimization process')
        parser.add_argument('--hic', type = str, nargs='+', help = 'Iteratively corrected Hi-C data')
        parser.add_argument('--chipseq', type = str, nargs = '+', help = 'ChIP-seq data')
        parser.add_argument('--np', type = int, default = 1, help = 'Number of processors')
        parser.add_argument('--resolution', type = int, default = 1, help = 'Resolution')
        parser.add_argument('--stepsize', type = float, default = 0.5, help = 'Step size to increment gamma parameter')
        parser.add_argument('--gamma_max', type = float, default = 4, help = 'Max gamma parameter')
        parser.add_argument('--hic_format', type = str, default = 'txt.gz', help = 'Hi-C matrices input format for armatus')
        parser.add_argument('--empty_row_imputation', action='store_true', help = 'Missing rows (and columns) imputation')
        parser.add_argument('--truncation',  action='store_true', help = 'Value truncation of Hi-C-matrix')
        parser.add_argument('--log2_hic', action='store_true', help = 'log2 transformation of Hi-C matrix')
        parser.add_argument('--log2_chip', action='store_true', help = 'log2 transformation of ChIP-Seq track')
        parser.add_argument('--zscore_chip', action='store_true', help = 'Z-score transformation of ChIP-Seq track')
        parser.set_defaults(empty_row_imputation=False)
        parser.set_defaults(truncation=False)
        parser.set_defaults(log2_hic=False)
        parser.set_defaults(log2_chip=False)
        parser.set_defaults(zscore_chip=False)
        args = parser.parse_args(sys.argv[2:])
        calculate.main(args, self.cfg, self.log)
        
        cpu_time = round(time.time()-start_time, 2)
        self.log.info('Execution time: {} sec'.format(cpu_time))
    
    def visualize(self):
        start_time = time.time()
        parser = argparse.ArgumentParser(description='Visualize results')
        parser.add_argument('--hic', type = str, help = 'Iteratively corrected Hi-C data')
        parser.add_argument('--tad', type = str, help = 'Set of Topologically associated domains')
        parser.add_argument('--region', type = str, help = 'Genome region')
        parser.add_argument('--resolution', type = int, help = 'Resolution')
        parser.add_argument('--chipseq', type = str, help = 'ChIP-seq data')
        parser.add_argument('--log2_chip', action='store_true', help = 'log2 transformation of ChIP-Seq track')
        parser.add_argument('--zscore_chip', action='store_true', help = 'Z-score transformation of ChIP-Seq track')
        args = parser.parse_args(sys.argv[2:])
        plotting.main(args, self.log)
        
        cpu_time = round(time.time()-start_time, 2)
        self.log.info('Execution time: {} sec'.format(cpu_time))


    
if __name__ == '__main__':
    optimalTAD()
