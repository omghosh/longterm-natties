#!/usr/bin/env bash
#SBATCH -J merge_PEAR
#SBATCH -p hns,dpetrov,normal,owners
#SBATCH -t 2-00:00
#SBATCH --array=1-95
#SBATCH --requeue
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -o SlurmFiles/slurm-%A_%a_%x.out
#SBATCH --mail-user=omghosh@stanford.edu
#SBATCH --mail-type=END

samples=$(sed -n "${SLURM_ARRAY_TASK_ID}"p samples_178.inp)


/home/groups/dpetrov/SOFTWARE/pear_0.9.11/bin/pear -f /scratch/groups/dpetrov/18032_178_yeast/${samples}_R1.fastq.gz -r /scratch/groups/dpetrov/18032_178_yeast/${samples}_R2.fastq.gz  -o merged_fastqs/${samples} -j $(nproc)



