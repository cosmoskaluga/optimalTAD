#!/bin/bash

cd armatus

if [ "$(uname)" == "Darwin" ]; then
    path_to_boost='/usr/local/Cellar/boost'  # assuming that boost was installed via brew. If not, change this to the actual path to boost directory.      
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    path_to_boost=$(whereis boost | awk '{print $2}')
elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ]; then
    echo "32 bits Windows NT platform is not supported"
elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW64_NT" ]; then
    echo "64 bits Windows NT platform is not supported"
fi

mpicxx -std=c++11 -w -L ${path_to_boost}/lib/ -I include/ -I ${path_to_boost}/include/ -O3 -o binaries/armatus src/*.cpp -lboost_iostreams -lboost_program_options -lboost_system
cp binaries/armatus ../optimalTAD/
cd ..
pip install .
