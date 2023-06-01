#!/bin/sh
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --job-name="adf"
#SBATCH --output=adf

module purge

cd $HOME

module load anaconda

source activate mlpot

cd mlpot-proj\src\python

python structure_analysis.py

source deactivate


