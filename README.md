[![Documentation Status](https://readthedocs.org/projects/optimaltad/badge/?version=latest)](https://optimaltad.readthedocs.io/en/latest)
[![DOI:10.1101/2023.03.06.531254](http://img.shields.io/badge/DOI-10.1101/2023.03.06.531254-B31B1B.svg)](https://doi.org/10.1101/2023.03.06.531254)

<a href="https://github.com/cosmoskaluga/optimalTAD"><img width="70%" src="https://github.com/cosmoskaluga/optimalTAD/blob/master/docs/images/optimalTAD_logo_square.png"></a>

This repo contains the source code of the algorithm for finding the optimal set of Topologically Associating Domains (TADs) using a combination of Hi-C and ChIP-seq/DNAme data. The algorithm is implemented in Python and is suitable for the identification of optimized TAD partitions across different resolutions in Drosophila and mammalian species.

## Getting Started

### Dependencies
To successfully run `optimalTAD` you will need to install the following dependencies:
- C++11
- Python 3
- MPICH2 or Open MPI (for parallel computing)
- boost (for macOS)

Also, one can set up a conda environment with all required dependencies installed using the following .yml files:

Linux users:
``` bash
conda env create -f environment_optimaltad_linux64.yml
```

macOS users:
``` bash
brew install boost # if boost libraries are not installed yet 
conda env create -f environment_optimaltad_osx.yml
``` 
Please note that the algorithm was tested on Linux and macOS operating systems only, therefore we can't guarantee that it will work on Windows as well.

### Installation
Clone this repo and install `optimalTAD` using pip:
``` bash
git clone https://github.com/cosmoskaluga/optimalTAD
cd optimalTAD
./install.sh
```

### Usage
To launch the algorithm type the following at the command line:
```bash
optimalTAD [-h] [--hic HIC [HIC ...]] [--chipseq CHIPSEQ [CHIPSEQ ...]] [--output OUTPUT] [--np NP] [--resolution RESOLUTION] [--stepsize STEPSIZE] [--gamma_max GAMMA_MAX] [--hic_format HIC_FORMAT] [--empty_row_imputation]
                  [--truncation] [--log2_hic] [--log2_chip] [--zscore_chip] [--balance | --no-balance] [--mammal] [--window_size_min WINDOW_SIZE_MIN] [--window_size_max WINDOW_SIZE_MAX]
```


Required and optional arguments:

    -h, --help                              Help message
    --hic HIC [HIC ...]                     Iteratively corrected Hi-C matrices in .hdf5 or .cool format
    --chipseq CHIPSEQ [CHIPSEQ ...]         Epigenetic data (ChIP-seq in .bedgraph or .bw format)
    --np [NP]                               Number of processors (=1)
    --resolution [RESOLUTION]               Resolution of Hi-C matrices (=1)
    --stepsize [STEPSIZE]                   Step size to increment gamma parameter in Armatus (=0.05)
    --gamma_max [GAMMA_MAX]                 Max value of the gamma parameter (=4)
    --hic_format [HIC_FORMAT]               Hi-C matrices input format for armatus (=txt.gz)
    --empty_row_imputation                  Empty line imputation (=False)
    --truncation                            Truncation of a Hi-C-matrix (=False)
    --log2_hic                              log2 transformation of Hi-C matrix (=False)
    --log2_chip                             log2 transformation of ChIP-seq values (=False)
    --zscore_chip                           Z-score transformation of ChIP-seq values(=False)
    --balance, --no-balance                 Hi-C matrix is iteratively normalized (--balance)
    --mammal                                Input data is derived from mammalian species (=False)
    --window_size_min [WINDOW_SIZE_MIN]     Minimal window size in insulation score method (for mammals only!)
    --window_size_max [WINDOW_SIZE_MAX]     Maximal window size in insulation score method (for mammals only!)  
    

All listed arguments can also be specified in the `config.ini` configuration file.

Both Hi-C and ChIP-Seq data are required for `optimalTAD` running. We strongly recommend you perform iterative correction ([Imakaev et al, 2012](https://www.nature.com/articles/nmeth.2148)) on your Hi-C data before running `optimalTAD`. ChIP-seq coverage track should be normalized by input and stored in .bedgraph or .bw file. No further preparation of ChIP-seq data is required, the algorithm will binarize coverage with respect to a chosen resolution of the Hi-C map and provide log2 and z-score transformation if needed.

`Note`: `optimalTAD` utilizes the two well-known TAD calling algorithms to produce all possible TAD sets and the choice of method depends on the species type input data originates from. The default method is [Armatus](https://github.com/kingsfordgroup/armatus), which is recommended for the analysis of Drosophila's topological domains. If you work with mammalian contact maps (`--mammals` argument), `optimalTAD` switches to the Insulation Score (IS) method from [cooltools](https://github.com/open2c/cooltools) package. 


### Running on test data
First, execute `test_data.sh` script:
```bash
chmod a+x ./test_data.sh
./test_data.sh
```
It will create `testdata` folder containing Hi-C and ChIP-seq data of Drosophila chromosome 2L. Next, run `optimalTAD`:

```bash
optimalTAD run
```
    
### Visualizing results
Hi-C data with the obtained optimal TAD set can be visualized using the function below:
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
    
 
## Documentation
`optimalTAD` documentation is available on [readthedocs](https://optimaltad.readthedocs.io/).

## How to cite
"optimalTAD: annotation of topologically associating domains based on chromatin marks enrichment." 
Dmitrii N. Smirnov, Anna D. Kononkova, Debra Toiber, Mikhail S. Gelfand, Ekaterina E. Khrameeva
bioRxiv 2023.03.06.531254; doi: https://doi.org/10.1101/2023.03.06.531254

## Manuscript
Scripts reproducing the main analyses from the manuscript can be found here:
https://github.com/cosmoskaluga/optimaltad_manuscript_analysis
