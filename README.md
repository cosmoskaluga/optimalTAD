[![Documentation Status](https://readthedocs.org/projects/optimaltad/badge/?version=latest)](https://optimaltad.readthedocs.io/en/latest)
[![BCH compliance](https://bettercodehub.com/edge/badge/cosmoskaluga/optimalTAD?branch=master)](https://bettercodehub.com/)
# optimalTAD
<a href="https://github.com/cosmoskaluga/optimalTAD"><img width="25%" src="https://github.com/cosmoskaluga/optimalTAD/blob/master/docs/optimalTAD_logo.png"></a>

Algorithm for finding the optimal set of Topologically Associating Domains based on the [Armatus](https://github.com/kingsfordgroup/armatus) software.

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
optimalTAD run [-h] [--hic HIC [HIC ...]] [--chipseq CHIPSEQ [CHIPSEQ ...]] [--np NP] [--resolution RESOLUTION] [--stepsize STEPSIZE] [--gamma_max GAMMA_MAX]
                   [--hic_format HIC_FORMAT] [--empty_row_imputation] [--truncation] [--log2_hic] [--log2_chip] [--zscore_chip]
```

Required and optional arguments:

    -h, --help                          Help message
    --hic HIC [HIC ...]                 Iteratively corrected Hi-C matrices in .hdf5 or .cool format
    --chipseq CHIPSEQ [CHIPSEQ ...]     Epigenetic data (ChIP-seq, ChIP-chip)
    --np [NP]                           Number of processors (=1)
    --resolution [RESOLUTION]           Resolution of Hi-C matrices (=1)
    --stepsize [STEPSIZE]               Step size to increment gamma parameter
    --gamma_max [GAMMA_MAX]             Max gamma parameter (=4)
    --hic_format [HIC_FORMAT]           Hi-C matrices input format for armatus (=txt.gz)
    --empty_row_imputation              Empty line imputation
    --truncation                        Truncation of a Hi-C-matrix (=False)
    --log2_hic                          log2 transformation of Hi-C matrix (=False)
    --log2_chip                         log2 transformation of epigenetic data (=False)
    --zscore_chip                       Z-score transformation of (=False)  
    
    
Hi-C data with the obtained optimal TAD set might be visualized using the function below:
```bash
optimalTAD visualize [-h] [--samplename SAMPLENAME] [--region REGION] [--resolution RESOLUTION] [--chipseq CHIPSEQ] [--log2_chip] [--zscore_chip] [--rnaseq RNASEQ]
```

with the following arguments:

    -h, --help                          Help message
    --samplename SAMPLENAME             Samplename of Hi-C data (for example, LacZ_1)
    --region REGION                     Genome region to plot (for example, chr2L:1,000,000-5,000,000)
    --resolution RESOLUTION             Resolution of Hi-C matrix (=1)
    --chipseq CHIPSEQ                   Path to epigenetic data
    --log2_chip                         log2 transformation of epigenetic data
    --zscore_chip                       Z-score transformation of epigenetic data
    --rnaseq RNASEQ                     Add additional track to the plot
    
 
### Documentation
The documentation for `optimalTAD` will be available on readthedocs.
                        

