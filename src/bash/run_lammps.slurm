#!/bin/sh
#SBATCH --time=8:00:00
#SBATCH --nodes=3
#SBATCH --mem=4000
#SBATCH --job-name="alo_test2"
#SBATCH --output=alo_test2_output.log

cd
cp -R $HOME/mlpot-proj/src/lammps_script/a_al2o3_input $WRKDIR
cd $WRKDIR

module purge

wget https://download.lammps.org/tars/lammps-23Jun2022.tar.gz
tar xf lammps-23Jun2022.tar.gz
cd lammps-23Jun2022

mkdir build
cd build
module load cmake gcc/8.4.0 openmpi/4.0.5

#build lammps:
cmake ../cmake -D PKG_KSPACE=yes -D PKG_MPIIO=yes -D CMAKE_INSTALL_PREFIX=../../lammps-mpi-23Jun2022
cmake --build .
make install

cd ../..

export PATH=$PATH:$PWD/lammps-mpi-23Jun2022/bin
cd $WRKDIR/a_al2o3_input

srun lmp -in a_al2o3_0_input.lammps

mv a_al2o3_0* $HOME

cd $WRKDIR
rm -r *

#for f in ./*input.lammps
#	do 
#		lmp -in f
#	done 

#compress .lammpstrj files with tar
#transfer to $HOME
