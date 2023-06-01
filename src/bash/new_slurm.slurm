#!/bin/sh
#SBATCH --time=8:00:00
#SBATCH --nodes=3
#SBATCH --mem=4000
#SBATCH --job-name="alo_test2"
#SBATCH --output=alo_test2_output.log

cd $WRKDIR

module purge
module load gcc openmpi

cd $WRKDIR/a_al2o3_input
export PATH=$PATH:/home/hiorthk1/lammps/build
srun lmp -in a_al2o3_0_input.lammps

#for f in ./*input.lammps
#	do 
#		lmp -in f
#	done 

#compress .lammpstrj files with tar
#transfer to $HOME
