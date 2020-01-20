#!/bin/bash
time=$(date +%s)

GREEN='\033[0;32m'
WHITE='\033[1;37m'
NC='\033[0m'
echo -e "${GREEN}TAD optimal size caller (2019)${NC}"
echo " "


#0. Cleaning old files
echo -e "${WHITE}0. Remove old files${NC}"
#rm -r Output/HiCmaps/*
rm -r Output/TADs/*
rm -r Output/NumeratedBorders/*
rm -r Output/Score/*
rm -r Output/Pictures/*
echo "     Done!"
echo " "

#1. Matrix extraction
echo -e "${WHITE}1. Start to prosess raw HiC map${NC}"
#python3.5 PyScripts/matrices_extraction.py True True
echo "     Done!"
echo " "

#2. Armatus calling
echo -e "${WHITE}2. TADs calling${NC}"
for fn in Output/HiCmaps/*
do
    echo "     ${fn:15}"
    mkdir Output/TADs/${fn:15}
    for name in Output/HiCmaps/${fn:15}/*
    do
        mydir=$(basename $name)

        if [ "$3" == "armatus" ]; then
            mpiexec -np $1 ./ArmatusParallel/src/armatus -r $2 -i $name -g 4.0 -o Output/TADs/${fn:15}/${mydir%%.*} -m -s 0.2
        else
            python3.5 PyScripts/lavaburst_calling.py $name Output/TADs/${fn:15}/${mydir%%.*} $2
        fi

    done
done
echo "     Done!"
echo " "



#3. TADs numeration
echo -e "${WHITE}3. TADs numeration${NC}"
for fn in Output/TADs/*
do
    echo "     ${fn:12}"
    python3.5 PyScripts/numeration.py $fn
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
    python3.5 PyScripts/stair_calling.py ${samples[$item]} ${repet[$item]}
done

echo "     Done!"
echo " "

#5. Best gamma Identification and Visualization
echo -e "${WHITE}5. Best gamma identification ${NC}"
python3.5 PyScripts/visualization.py
echo "     Done!"
echo " "

echo "CPU time: $(($(date +%s)-$time))s"
