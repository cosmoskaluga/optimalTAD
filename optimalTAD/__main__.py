import argparse
import logging
import sys
import time
from . import calculate
from . import logger


class optimalTAD:
    def __init__(self):
        self.log = logger.initialize_logger()
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
        parser.add_argument('--resolution', type = int, default = 1, help = 'Resolution of Hi-C matrices')
        parser.add_argument('--stepsize', type = float, default = 0.5, help = 'Step size to increment gamma parameter')
        parser.add_argument('--gamma_max', type = float, default = 4, help = 'Max gamma parameter')
        parser.add_argument('--hic_format', type = str, default = 'txt.gz', help = 'Hi-C matrices input format for armatus')
        args = parser.parse_args(sys.argv[2:])
        calculate.main(args, self.log)
        self.log.info('Execution time: {} sec'.format(time.time()-start_time))
    
    def visualize(self):
        parser = argparse.ArgumentParser(description='Visualize results')
        args = parser.parse_args(sys.argv[2:])



    
if __name__ == '__main__':
    optimalTAD()
