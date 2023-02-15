#!/bin/sh
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --mem=2000
#SBATCH --job-name="lammps-test"
#SBATCH --output=test1_output.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hiorthk1@aalto.fi


cp -R mlpot-proj/src/lammps_script/a_al2o3_input $WRKDIR
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
make -j 2
make install


cd ../..

#add installed lammps to the executable search path
export PATH = $PATH:$PWD/lammps-mpi-23Jun2022/bin

cd $WRKDIR/a_al2o3_input

lmp -in a_al2o3_0_input.lammps

mv *.lammpstrj $HOME
mv *.log $HOME


#for f in ./*input.lammps
#	do 
#		lmp -in f
#	done 

#compress .lammpstrj files with tar
#transfer to $HOME
