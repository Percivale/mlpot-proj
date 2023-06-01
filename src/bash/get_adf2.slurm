#!/bin/sh
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --mem=1000
#SBATCH --job-name="adf"
#SBATCH	--output=adf.log

module purge

module load anaconda
source activate mlpot

cd $WRKDIR
cd mlpot/python

python structure_analysis.py

conda deactivate
