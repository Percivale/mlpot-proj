#!/bin/sh
#SBATCH --partition=short
#SBATCH --account=share-nv-fys
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1000
#SBATCH --job-name="git_test"

WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}

module purge
module load git/2.36.0-GCCcore-11.3.0-nodocs


#git clone "https://github.com/Percivale/mlpot-proj"
cd mlpot-proj
git pull origin main

sed 's/wrld/world' doc/test.txt

git add .
git commit -m "Idun test"
git push origin main


