#!/usr/bin/env bash
#SBATCH -J count_bcs
#SBATCH -p hns,dpetrov,normal
#SBATCH -t 0-02:00
#SBATCH --array=1-385
#SBATCH --requeue
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH -o SlurmFiles/slurm-%A_%a_%x.out
#SBATCH --mail-user=omghosh@stanford.edu
#SBATCH --mail-type=END


. /home/users/omghosh/.bashrc
conda activate natty_env

samples=$(sed -n "${SLURM_ARRAY_TASK_ID}"p samples_bc_counting.inp)

python3 optimized_bc_counting_Jan2025.py $samples