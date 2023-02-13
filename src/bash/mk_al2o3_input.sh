#!/bin/bash
for i in {0..10}
do
	echo "initializing a_al2o3_${i}_input.lammps"
	cp a_al2o3_input.lammps a_al2o3_input/a_al2o3_${i}_input.lammps
	sed -i "s/NUMBER_ID/${i}/" a_al2o3_input/a_al2o3_${i}_input.lammps
	sed -i "s/SEED1/$(($i+1))/" a_al2o3_input/a_al2o3_${i}_input.lammps
	sed -i "s/SEED2/$(($i+2))/" a_al2o3_input/a_al2o3_${i}_input.lammps
done
