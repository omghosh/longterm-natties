#!/usr/bin/env bash
#SBATCH -J process_reads_167
#SBATCH -p hns,dpetrov,normal
#SBATCH -t 2-00:00
#SBATCH --array=1-325%100
#SBATCH --requeue
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH -o SlurmFiles/slurm-%A_%a_%x.out
#SBATCH --mail-user=omghosh@stanford.edu
#SBATCH --mail-type=END


. /home/users/omghosh/.bashrc
conda activate natty_env

samples=$(sed -n "${SLURM_ARRAY_TASK_ID}"p samples_167.inp)

python3 process_merged_reads.py $samples
