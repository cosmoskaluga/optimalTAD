# Armatus Parallel
A parallel version of the Armatus software (https://github.com/kingsfordgroup/armatus)

## MPICH2 installation
Install MPICH2 (Ubuntu / Linux Mint / Debian)
```bash
sudo apt update
sudo apt install mpich2
```
Install MPICH2 on Mac OSX
```bash
brew install mpich2
```
## Building
To compile the code type the following commands in the terminal:
```bash
cmake CMakeLists.txt
make
```
## Run
```bash
./run_armatus.sh
```
