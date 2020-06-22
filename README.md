# optimalTAD
This repository is by Dmitrii Smirnov and Ekaterina Khrameeva and contains the source code of the algorithm of finding the optimal set of Topologically Associating Domains via [Armatus](https://github.com/kingsfordgroup/armatus) software.

## Getting Started

### Dependencies
- Linux or macOS
- Python 3
- MPICH2 

### Installation
Clone this repo:
``` bash
git clone https://github.com/cosmoskaluga/optimalTAD
cd optimalTAD
pip install .
```
### Usage
```bash
python3 -m optimalTAD [--hic HIC] [--chipseq CHIPSEQ] [--np NP] [--resolution RESOLUTION] [--stepsize STEPSIZE] [--gamma_max GAMMA_MAX] [--hic_format HIC_FORMAT]
```
Optional arguments:

    -h, --help            show this help message and exit
    --hic HIC
    --chipseq CHIPSEQ
    --np NP
    --resolution RESOLUTION
    --stepsize STEPSIZE
    --gamma_max GAMMA_MAX
    --hic_format HIC_FORMAT

