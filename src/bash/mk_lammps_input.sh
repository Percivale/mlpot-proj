#!/bin/bash
#run from dir with input files

for ((i = 1; i<=10;i++))
do
	echo "initializing a_al2o3${i}_input.lammps"
	
	cp a_al2o3_input.lammps a_al203_input/a_al2o3_${i}_input.lammps
	
	

done