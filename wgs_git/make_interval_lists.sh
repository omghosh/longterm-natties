#!/usr/bin/env bash
#SBATCH -J gmake_intervals_list 
#SBATCH -p hns,dpetrov,normal,owners
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -t 2-00:00
#SBATCH --array=1-11
#SBATCH --mem-per-cpu=2G
#SBATCH --requeue
#SBATCH -o SlurmFiles/slurm-%A_%a_%x.out

# Set the reference strain variable dynamically based on the SLURM_ARRAY_TASK_ID
ref_strain=$(sed -n "$SLURM_ARRAY_TASK_ID"p all_reference_strains.inp)

# Base directory for reference genomes
ref_base_dir="/scratch/groups/dpetrov/wgs_natties/reference_assemblies/GENOMES_ASSEMBLED"

# Construct the reference genome path dynamically
ref_genome="$ref_base_dir/${ref_strain}.re.fasta"

# Output file for intervals
output_file="intervals_${ref_strain}.list"

# Extract lines that start with '>' from the reference genome file and write to the output file
grep "^>" "$ref_genome" | sed 's/^>//' > "$output_file"

echo "Generated intervals file: $output_file"