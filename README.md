[![Documentation Status](https://readthedocs.org/projects/optimaltad/badge/?version=latest)](https://optimaltad.readthedocs.io/en/latest)
[![BCH compliance](https://bettercodehub.com/edge/badge/cosmoskaluga/optimalTAD?branch=master)](https://bettercodehub.com/)

<a href="https://github.com/cosmoskaluga/optimalTAD"><img width="70%" src="https://github.com/cosmoskaluga/optimalTAD/blob/master/docs/optimalTAD_logo_square.png"></a>

This repo contains source code of the algorithm for finding the optimal set of Topologically Associating Domains (TADs) using a combination of Hi-C and ChIP-seq data. The algorithm is implemented in Python and incrporates [Armatus](https://github.com/kingsfordgroup/armatus) tool for calling TADs across different resolutions.

## Getting Started

### Dependencies
To succesfully run `optimalTAD` you will need to install following dependencies:
- C++11
- Python 3
- MPICH2 or Open MPI (for parallel computing)
- boost (for macOS)

Please note that the algorithm was tested on Linux and macOS operating systems only, therefore we can't guarantee that it will work on Windows as well.

### Installation
Clone this repo and install `optimalTAD` using pip:
``` bash
git clone https://github.com/cosmoskaluga/optimalTAD
cd optimalTAD
pip install .
```

### Usage
To launch the algorithm type the following at the command line:
```bash
python3 -m optimalTAD run [-h] [--hic HIC [HIC ...]] [--chipseq CHIPSEQ [CHIPSEQ ...]] [--np NP] [--resolution RESOLUTION] [--stepsize STEPSIZE] [--gamma_max GAMMA_MAX] [--hic_format HIC_FORMAT] [--empty_row_imputation] [--truncation] [--log2_hic] [--log2_chip] [--zscore_chip]
```



Required and optional arguments:

    -h, --help                          Help message
    --hic HIC [HIC ...]                 Iteratively corrected Hi-C matrices in .hdf5 or .cool format
    --chipseq CHIPSEQ [CHIPSEQ ...]     Epigenetic data (ChIP-seq in .bedgraph or .bw format)
    --np [NP]                           Number of processors (=1)
    --resolution [RESOLUTION]           Resolution of Hi-C matrices (=1)
    --stepsize [STEPSIZE]               Step size to increment gamma parameter in Armatus (=0.05)
    --gamma_max [GAMMA_MAX]             Max value of the gamma parameter (=4)
    --hic_format [HIC_FORMAT]           Hi-C matrices input format for armatus (=txt.gz)
    --empty_row_imputation              Empty line imputation (=False)
    --truncation                        Truncation of a Hi-C-matrix (=False)
    --log2_hic                          log2 transformation of Hi-C matrix (=False)
    --log2_chip                         log2 transformation of ChIP-seq values (=False)
    --zscore_chip                       Z-score transformation of ChIP-seq values(=False)  
    

All listed arguments can also be specified in the `config.ini` configuration file.

Both Hi-C and ChIP-Seq data are required for `optimalTAD` running. We strongly recommend you to perform iterative correction ([Imakaev et al, 2012](https://www.nature.com/articles/nmeth.2148)) on your Hi-C data before running `optimalTAD`. ChIP-seq coverage track should be normalized by input and stored in .bedgraph or .bw file. No further preparation of ChIP-seq data is required, algorithm will binarize coverage with respect to chosen resolution of Hi-C map and provide log2 and z-score transformation if needed.
    
### Visualizing results
Hi-C data with the obtained optimal TAD set can be visualized using the function below:
```bash
python3 -m optimalTAD visualize [-h] [--samplename SAMPLENAME] [--region REGION] [--resolution RESOLUTION] [--chipseq CHIPSEQ] [--log2_chip] [--zscore_chip] [--rnaseq RNASEQ]
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
    
 
## Documentation
The documentation for `optimalTAD` will be available on readthedocs.
                        

