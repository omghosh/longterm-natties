#!/usr/bin/env bash
#SBATCH -J gtak_mergeall
#SBATCH -p hns,dpetrov,normal,owners
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -t 2-00:00
#SBATCH --array=1-11
#SBATCH --mem-per-cpu=2G
#SBATCH --requeue
#SBATCH -o SlurmFiles/slurm-%A_%a_%x.out

ml system
ml ncurses/6.0
ml biology
ml bwa/0.7.17
ml java

# Set the reference strain variable dynamically based on the SLURM_ARRAY_TASK_ID
ref_strain=$(sed -n "$SLURM_ARRAY_TASK_ID"p all_reference_strains.inp)

# Base directory for reference genomes
ref_base_dir="/scratch/groups/dpetrov/wgs_natties/reference_assemblies/GENOMES_ASSEMBLED"

# Construct the reference genome path dynamically
ref_genome="$ref_base_dir/${ref_strain}.re.fasta"

# Construct the sample map path dynamically
sample_map="yeast_samples_${ref_strain}.sample_map"

# Point to data directory and output directory
dir="RawData/"
out_dir="Output/"

# Construct the output VCF path dynamically
output_vcf="${out_dir}WGS_AllSamples_${ref_strain}.vcf.gz"

echo "Using reference genome: $ref_genome"
echo "Using sample map: $sample_map"
echo "Output VCF will be: $output_vcf"

### Merge GVCF files

# (version for many files) - need to make a sample map
/home/groups/dpetrov/SOFTWARE/gatk-4.2.0.0/gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport -L intervals_${ref_strain}.list --sample-name-map "$sample_map" --genomicsdb-workspace-path "${out_dir}genomicsdb_${ref_strain}"

### Call genotypes for all samples at the same time
/home/groups/dpetrov/SOFTWARE/gatk-4.2.0.0/gatk --java-options "-Xmx4g" GenotypeGVCFs -R "$ref_genome" -V gendb://"$out_dir"/genomicsdb_${ref_strain} -O "$output_vcf"