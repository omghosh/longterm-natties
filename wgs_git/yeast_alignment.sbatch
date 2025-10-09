#!/usr/bin/env bash
#SBATCH -J bwa
#SBATCH -p hns,dpetrov,normal,owners
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -t 2-00:00
#SBATCH --array=1-45,71-139,163-182
#SBATCH --mem-per-cpu=8G
#SBATCH --requeue
#SBATCH -o SlurmFiles/slurm-%A_%a_%x.out

bn=$(sed -n "$SLURM_ARRAY_TASK_ID"p wgs_merged_samples.inp)


ml system
ml ncurses/6.0
ml biology
ml bwa/0.7.17
# Try to load a newer Java version first, fall back to default if not available
ml java/11 2>/dev/null || ml java/1.8.0_latest 2>/dev/null || ml java
ml samtools

# Create associative array from external file
declare -A ref_genomes
while IFS= read -r line; do
    eval "$line"
done < associative_array.sh

# Print the entire associative array for debugging
echo "Associative array contents:"
for key in "${!ref_genomes[@]}"; do
    echo "$key => ${ref_genomes[$key]}"
done



ref_genome=${ref_genomes[$bn]}

if [ -z "$ref_genome" ]; then
    echo "Error: Reference genome not defined for sample $bn"
    exit 1
fi

dir="/scratch/groups/dpetrov/wgs_natties/merged_fastq_files/"
out_dir="Output/"

# Create output directory if it doesn't exist
mkdir -p "$out_dir"

# Set up a temporary directory for large operations
export TMPDIR="/scratch/groups/dpetrov/tmp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"

# Clean up temp directory on exit
trap "rm -rf $TMPDIR" EXIT

### Index checks (your existing code is good)
if [ ! -f "${ref_genome}.bwt" ] || [ ! -f "${ref_genome}.sa" ]; then
    echo "Creating BWA index..."
    if ! bwa index "$ref_genome"; then
        echo "Error: BWA indexing failed"
        exit 1
    fi
    echo "BWA indexing completed successfully"
else
    echo "BWA index already exists, skipping indexing step"
fi

if [ ! -f "${ref_genome%.*}.dict" ]; then
    echo "Creating GATK dictionary..."
    /home/groups/dpetrov/SOFTWARE/gatk-4.2.0.0/gatk CreateSequenceDictionary -R "$ref_genome"
else
    echo "GATK dictionary already exists, skipping"
fi

if [ ! -f "${ref_genome}.fai" ]; then
    echo "Creating samtools index..."
    samtools faidx "$ref_genome"
else
    echo "Samtools index already exists, skipping"
fi

### Align to reference
bwa mem -t 16 "$ref_genome" "$dir$bn"_R1.fastq.gz "$dir$bn"_R2.fastq.gz | samtools view -bS - | samtools sort > "$out_dir$bn".sorted.bam

### Fix read groups - ADD MEMORY AND GC FLAGS
java -Xmx16g -XX:+UseG1GC -XX:MaxGCPauseMillis=200 \
     -Djava.io.tmpdir="$TMPDIR" \
     -jar /home/groups/dpetrov/SOFTWARE/picard-tools/picard.jar AddOrReplaceReadGroups \
     I="$out_dir$bn".sorted.bam \
     O="$out_dir$bn".sorted.rg.bam \
     RGID="$bn" \
     RGLB=lib"$bn" \
     RGPL=illumina \
     RGPU=unit"$bn" \
     RGSM="$bn" \
     TMP_DIR="$TMPDIR"

samtools sort "$out_dir$bn".sorted.rg.bam -o "$out_dir$bn".sorted.rg.coord.bam

samtools index "$out_dir$bn".sorted.rg.coord.bam  # Index the coord.bam, not the unsorted one

### Mark duplicates - THIS IS THE CRITICAL PART
java -Xmx32g -XX:+UseG1GC -XX:MaxGCPauseMillis=200 \
     -Djava.io.tmpdir="$TMPDIR" \
     -jar /home/groups/dpetrov/SOFTWARE/picard-tools/picard.jar MarkDuplicates \
     I="$out_dir$bn".sorted.rg.coord.bam \
     O="$out_dir$bn".sorted.rg.md.bam \
     M="$out_dir$bn".md.metrics.txt \
     MAX_RECORDS_IN_RAM=250000 \
     MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
     SORTING_COLLECTION_SIZE_RATIO=0.25 \
     TMP_DIR="$TMPDIR" \
     VALIDATION_STRINGENCY=LENIENT

# Check if MarkDuplicates succeeded before proceeding
if [ ! -f "$out_dir$bn".sorted.rg.md.bam ]; then
    echo "Error: MarkDuplicates failed, output file does not exist"
    exit 1
fi

samtools index "$out_dir$bn".sorted.rg.md.bam

### Call variants with GATK - INCREASE MEMORY HERE TOO
/home/groups/dpetrov/SOFTWARE/gatk-4.2.0.0/gatk --java-options "-Xmx16g -XX:+UseG1GC" \
     HaplotypeCaller \
     -R "$ref_genome" \
     -I "$out_dir$bn".sorted.rg.md.bam \
     -ERC GVCF \
     -O "$out_dir$bn".g.vcf.gz \
     --tmp-dir "$TMPDIR"

echo "Pipeline completed successfully for sample: $bn"
