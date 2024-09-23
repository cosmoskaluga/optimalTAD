#!/bin/bash

mkdir -p testdata/fly && cd $_

wget -O LacZ1_2L.bedgraph "https://www.dropbox.com/s/trws51e5wdgspn7/LacZ1_2L.bedgraph?dl=1"
wget -O LacZ1_2L.hdf5 "https://www.dropbox.com/s/bpjmqfalmgp8c5j/LacZ1_2L.hdf5?dl=1"

mkdir ../mammal && cd $_
wget -O mammal_chr1.bedgraph "https://www.dropbox.com/scl/fi/qek9gne57u8ki7479dq4h/mammal_chr1.bedgraph?rlkey=qj8d3c19u46ytz334ie47nrdl&dl=1"
wget -O mammal_chr1.cool "https://www.dropbox.com/scl/fi/urvqcmf2zf5t1m1lmdsjn/mammal_chr1.cool?rlkey=vdymiktwdvbc5sd17jvfzery4&dl=1"






