# optimalTAD
This repository is by Dmitrii Smirnov, Ekaterina Khrameeva and contains the source code of the algorithm of finding the optimal set of Topologically Associating Domains via [Armatus](https://github.com/kingsfordgroup/armatus) software.

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
```
### Usage
```bash
./run.sh np res caller
```
The command line arguments: 
- np: number of CPUs (for Armatus only),
- res: resolution,
- caller: armatus, lavaburst
