#!/bin/bash
time=$(date +%s)

GREEN='\033[0;32m'
WHITE='\033[1;37m'
NC='\033[0m'
echo -e "${GREEN}TAD optimal size caller (version 0.1)${NC}"
echo " "

usage() { echo "Usage: $0 [-c <armatus|matreshka>] [-p <cores>] [-r <resolution>] [-s <stepsize>]" 1>&2; exit 1; }

while getopts ":c:p:r:s:h" o; do
    case "${o}" in
        c)
            caller=${OPTARG}
            if [ $caller != "armatus" ] && [ $caller != "matreshka" ]; then
                usage
            fi
            ;;
        p)
            np=${OPTARG}
            ;;
        r)
            resolution=${OPTARG}
            ;;
        s)
            stepsize=${OPTARG}
            ;;
        h)
            echo "Usage:"
            echo "      -c | name of caller (armatus or matreshka)"
            echo "      -p | number of CPU"
            echo "      -r | resolution (armatus parameter)"
            echo "      -s | stepsize (armatus parameter)"
            exit 1;
            ;;
        *)
            echo "No reasonable options found!"
            usage
            ;;
    esac
done

#0. Cleaning old files
echo -e "${WHITE}0. Remove old files${NC}"
rm -rf Output/HiCmaps/* || true
rm -rf Output/TADs/* || true
rm -rf Output/NumeratedBorders/* || true
rm -rf Output/Score/* || true
rm -rf Output/Pictures/* || true
echo "     Done!"
echo " "


#1. Matrix extraction
echo -e "${WHITE}1. Start to prosess raw HiC map${NC}"
python3 PyScripts/matrices_extraction.py True True
echo "     Done!"
echo " "


#2. Armatus calling
echo -e "${WHITE}2. TAD calling${NC}"
for path_to_folder in Output/HiCmaps/*
do
    foldername=${path_to_folder##*/}
    echo "     ${foldername}"
    mkdir Output/TADs/${foldername}
    for name in Output/HiCmaps/${foldername}/*
    do
        output_dir=$(basename $name)

        if [ "$caller" = "armatus" ]; then
            gamma_max=4.0
            path_to_armatus="ArmatusParallel/src/armatus"
            path_to_output="Output/TADs/${foldername}/${output_dir%%.*}"

            mpiexec -np $np ./${path_to_armatus} -r $resolution \
                                                    -i $name \
                                                    -g $gamma_max \
                                                    -o $path_to_output \
                                                    -s $stepsize \
                                                    -m
        else
            python3 PyScripts/lavaburst_calling.py $name Output/TADs/${foldername}/${output_dir%%.*} ${resolution}
        fi
    done
done
echo "     Done!"
echo " "


#3. TADs numeration
echo -e "${WHITE}3. TADs numeration${NC}"
for filename in Output/TADs/*
do
    echo "     ${filename##*/}"
    python3 PyScripts/numeration.py $filename
done
echo "     Done!"
echo " "


#4. Stairs calling
echo -e "${WHITE}4. Stairs calling${NC}"

file="ChIP_Seq/Repet.txt"
samples=()
repet=()

while read line
do
    samples+=("$(cut -d"." -f1 <<<$line)")
    repet+=("$(cut -d" " -f2 <<<$line)")
done <"$file"

len=${#samples[*]}

for (( item = 0; item < len; item++ ))
do
    echo "     ${samples[$item]}"
    python3 PyScripts/stair_calling.py ${samples[$item]} ${repet[$item]}
done

echo "     Done!"
echo " "

#5. Best gamma Identification and Visualization
echo -e "${WHITE}5. Best gamma identification ${NC}"
python3 PyScripts/visualization.py
echo "     Done!"
echo " "

echo "CPU time: $(($(date +%s)-$time))s"
