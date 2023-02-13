#!/bin/sh
#SBATCH --partition=short
#SBATCH --account=share-nv-fys
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1G
#SBATCH --job-name="test1"
#SBATCH --output=test1_output.log

WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}

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


module load git/2.21.0
git clone "repo_name"

cd git_repo/
#do stuff

git push origin main #can I do this?
