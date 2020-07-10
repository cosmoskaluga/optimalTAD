ARMATUS=src/armatus
echo $gamma

mpiexec -np 4 ./$ARMATUS -r 20000 -i examples/chr2L_control.txt.gz -g 5.2 -o output/test_gamma${gamma} -m -s 0.5

