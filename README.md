[![BCH compliance](https://bettercodehub.com/edge/badge/cosmoskaluga/optimalTAD?branch=master)](https://bettercodehub.com/)


# optimalTAD
Algorithm for finding the optimal set of Topologically Associating Domains based on [Armatus](https://github.com/kingsfordgroup/armatus) software.

## Getting Started

### Dependencies
You will need the following dependencies:
- C++11
- Python 3
- MPICH2 (for parallel computing)
- boost (for macOS)

Please note that the algorithm was tested on Linux and macOS operating systems only, therefore we can't guarantee correctness of its execution on Windows.

### Installation
Clone this repo and install `optimalTAD` using pip:
``` bash
git clone https://github.com/cosmoskaluga/optimalTAD
cd optimalTAD
pip install .
```

### Input format
Both Hi-C and ChIP-Seq data are required for `optimalTAD` running. Hi-C data should be iteratively corrected ([Imakaev et al, 2012](https://www.nature.com/articles/nmeth.2148)) and contained in .hdf5 file. Bed format for acetylation data is needed. Please check the documentation for more details.

### Usage
To launch the algorithm just type the following at the command line:
```bash
optimalTAD [--hic HIC] [--chipseq CHIPSEQ] [--np NP] [--resolution RESOLUTION] [--stepsize STEPSIZE] [--gamma_max GAMMA_MAX] [--hic_format HIC_FORMAT]
```

Optional arguments:

    -h, --help                     Help message
    --hic [HIC]                    Iteratively corrected Hi-C matrices in hdf5 format
    --chipseq [CHIPSEQ]            ChIP-Seq data
    --np [NP]                      Number of processors (=1)
    --resolution [RESOLUTION]      Resolution of Hi-C matrices (=1)
    --stepsize [STEPSIZE]          Step size to increment gamma parameter
    --gamma_max [GAMMA_MAX]        Max gamma parameter (=4)
    --hic_format [HIC_FORMAT]      Hi-C matrices input format for armatus (=txt.gz)
 
### Documentation
The documentation for `optimalTAD` will be available on readthedocs.
                        

