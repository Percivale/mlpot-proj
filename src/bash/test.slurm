#!/bin/sh
#SBATCH --partition=CPUQ
#SBATCH --account=kajahh
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1000
#SBATCH --job-name="git_test"
#SBATCH --output=test-srun.out

WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}

echo "we are running from this directory: $SLURM_SUBMIT_DIR"
echo " the name of the job is: $SLURM_JOB_NAME"
echo "Th job ID is $SLURM_JOB_ID"
echo "The job was run on these nodes: $SLURM_JOB_NODELIST"
echo "Number of nodes: $SLURM_JOB_NUM_NODES"
echo "We are using $SLURM_CPUS_ON_NODE cores"
echo "We are using $SLURM_CPUS_ON_NODE cores per node"
echo "Total of $SLURM_NTASKS cores"

module purge
module load git/2.36.0-GCCcore-11.3.0-nodocs


git clone "https://github.com/Percivale/mlpot-proj"
cd mlpot-proj

sed "ReadMe.md"

