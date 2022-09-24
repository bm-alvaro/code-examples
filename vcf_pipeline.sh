#!/usr/bin/env bash

# Number of desired cpus:
#SBATCH --cpus-per-task=52

# Amount of RAM needed for this job:
#SBATCH --mem=180gb

# The time the job will be running:
#SBATCH --time=6-23:59:00

# Set output and error files
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

# Loading modules:
module load bwa/0.7.12
module load samtools/1.13
module load gatk/4.2.2.0
module load picard/2.7.1
module load bcftools/1.4

SCRIPT_NAME=$(basename $0)

RAW=/mnt/scratch/users/bio_293_ual/alvarobm/BGI_data/test_raw
RESOURCES=/mnt/scratch/users/bio_293_ual/alvarobm/BGI_data/new_run_10.1002_pld3.404/RESOURCES
OUTPUT=/mnt/scratch/users/bio_293_ual/alvarobm/BGI_data/new_run_10.1002_pld3.404/OUTPUT
PICARD=/mnt/home/soft/picard/programs/x86_64/2.7.1
THREADS=52
CORES=52

echo "########################################################################################"
echo "############################# Mapping and variant calling ##############################"  
echo "########################################################################################"
echo "                              Alvaro Benitez Mateo, 2022"

echo "The pathing to be use is:
$RAW as the folder containing raw data
$RESOURCES containing the reference genome
$OUTPUT as the directory to stream the data output"

# Creating output subdirectories

mkdir -p $OUTPUT/bwa/logs
mkdir -p $OUTPUT/samtools/logs
mkdir -p $OUTPUT/gatk-individual/logs
mkdir -p $OUTPUT/gatk-individual/samtools
mkdir -p $OUTPUT/gatk-joint/logs

# Defining directory variables

OUT_BWA=$OUTPUT/bwa
LOGS_BWA=$OUT_BWA/logs
OUT_SAMTOOLS=$OUTPUT/samtools
LOGS_SAMTOOLS=$OUT_SAMTOOLS/logs
OUT_GATK_IND=$OUTPUT/gatk-individual
LOGS_GATK_IND=$OUT_GATK_IND/logs
OUT_GATKSM_IND=$OUT_GATK_IND/samtools
OUT_GATK_JOINT=$OUTPUT/gatk-joint
LOGS_GATK_JOINT=$OUT_GATK_JOINT/logs

################################################

echo "BWA (indexing step) ---- STARTING"

echo "Indexing reference genome"

bwa index -p $RESOURCES/cpepo_genome_v4.1 $RESOURCES/Cpepo_genome_v4.1.fa

echo "Indexing done"

echo "BWA (indexing step) ---- FINISHED"

echo "Alignment with bwa"

echo "BWA (mapping step) ---- STARTING"

for R1 in $RAW/*1.fastq.gz

    do

        FILE_NAME=$(basename $R1 | sed "s/1.fastq.gz//" - )
        FILE_NAME_R1=$(basename $R1 | sed "s/.fastq.gz//" - )
        R2=$(echo $R1 | sed "s/1.fastq.gz/2.fastq.gz/" - )
        FILE_NAME_R2=$(basename $R2 | sed "s/_2.fastq.gz//" - )

        echo "##################################################"

        date
        echo "Left pair: $FILE_NAME_R1"
        echo "Right pair: $FILE_NAME_R2"

        echo "##################################################"

        SAMPLE=$FILE_NAME
        PLATFORM="ILLUMINA"
        LIBRARY="Lib-1"
        TPU="none"
        RGID=${SAMPLE}_${LIBRARY}
        RG="@RG\tID:${RGID}\tSM:${SAMPLE}\tPL:${PLATFORM}\tLB:${LIBRARY}\tPU:${TPU}"

        time bwa mem \
            -R $RG \
            $RESOURCES/cpepo_genome_v4.1 \
            $R1 $R2 \
            -t $THREADS > $OUT_BWA/$FILE_NAME.bwa_mem_align.sam 2>> $LOGS_BWA/bwa_mem.stderr

        echo "Processing of the file $FILE_NAME is finished :-) !!!"
        echo "##################################################"

done

echo "BWA (mapping step) ---- FINISHED"

# Second part: SAMtools based variant-calling

echo "SAMtools (reference indexing step) ---- STARTING"

samtools faidx $RESOURCES/Cpepo_genome_v4.1.fa > $RESOURCES/Cpepo_genome_v4.1.fa.fai

echo "SAMtools (reference indexing step) ---- FINISHED"

echo "SAMtools (sorting step) ---- STARTING"

for SAM in $OUT_BWA/*.sam

    do

        FILE_NAME=$(basename $SAM | sed "s/.bwa_mem_align.sam//" - )

        echo "##################################################"

        date
        echo "Sample: $FILE_NAME"

        echo "##################################################"

        time samtools sort -@ $THREADS -m 2G -O bam -T $OUT_SAMTOOLS/$FILE_NAME.tmp -o $OUT_SAMTOOLS/$FILE_NAME.sorted.bam $SAM 2>> $LOGS_SAMTOOLS/sort.stderr

        echo "Processing of the file $FILE_NAME is finished :-) !!!"
        echo "##################################################"

done

echo "SAMtools (sorting step) ---- FINISHED"

echo "PICARD (identifying duplicates step) ---- STARTING"

for BAM in $OUT_SAMTOOLS/*.sorted.bam

    do
        
        FILE_NAME=$(basename $BAM | sed "s/.sorted.bam//" - )

        echo "##################################################"

        date
        echo "Sample: $FILE_NAME"

        echo "##################################################"

        time java -jar $PICARD/picard.jar MarkDuplicates \
            I=$BAM \
            O=$OUT_SAMTOOLS/$FILE_NAME.marked_duplicates.bam \
            M=$OUT_SAMTOOLS/$FILE_NAME.marked_dup_metrics.txt 2>> $LOGS_SAMTOOLS/MarkDuplicates.stderr

        echo "Processing of the file $FILE_NAME is finished :-) !!!"
        echo "##################################################"

done

echo "PICARD (identifying duplicates step) ---- FINISHED"

echo "SAMtools (indexing bam-files step) ---- STARTING"

for BAM in $OUT_SAMTOOLS/*.marked_duplicates.bam

    do

        FILE_NAME=$(basename $BAM | sed "s/.marked_duplicates.bam//" - )

        echo "##################################################"

        date
        echo "Sample: $FILE_NAME"

        echo "##################################################"

        time samtools index -b -@ $THREADS $BAM 2>>$LOGS_SAMTOOLS/samtools_index.stderr

        echo "Processing of the file $FILE_NAME is finished :-) !!!"
        echo "##################################################"

done

echo "SAMtools (indexing bam-files step) ---- FINISHED"

echo "BCFtools (genotype likelihoods step) ---- STARTING"

for BAM in $OUT_SAMTOOLS/*.marked_duplicates.bam

    do

        FILE_NAME=$(basename $BAM | sed "s/.marked_duplicates.bam//" - )

        echo "##################################################"

        date
        echo "Sample: $FILE_NAME"

        echo "##################################################"

        time bcftools mpileup -Ou -f $RESOURCES/Cpepo_genome_v4.1.fa $BAM | bcftools call -vmO z > $OUT_SAMTOOLS/$FILE_NAME.vcf 2>>$LOGS_SAMTOOLS/bcftools_mpileup.stderr

        echo "Processing of the file $FILE_NAME is finished :-) !!!"
        echo "##################################################"

done

echo "BCFtools (genotype likelihoods step) ---- FINISHED"

# Third part: GATK based variant-calling [Individual variant-calling]

echo "GATK (reference indexing step) ---- STARTING"

java -jar $PICARD/picard.jar CreateSequenceDictionary \ 
      R=$RESOURCES/Cpepo_genome_v4.1.fa \ 
      O=$RESOURCES/Cpepo_genome_v4.1.dict

echo "GATK (reference indexing step) ---- FINISHED"

echo "GATK (variant calling step (default mode)) ---- STARTING"

for BAM in $OUT_SAMTOOLS/*.marked_duplicates.bam

    do

        FILE_NAME=$(basename $BAM | sed "s/.marked_duplicates.bam//" - )

        echo "##################################################"

        date
        echo "Sample: $FILE_NAME"

        echo "##################################################"

        time  gatk --java-options "-Xmx4g" HaplotypeCaller  \
                -R $RESOURCES/Cpepo_genome_v4.1.fa \
                -I $BAM \
                -O $OUT_GATK_IND/$FILE_NAME.g.vcf.gz > $LOGS_GATK_IND/$FILE_NAME.stdout  2> $LOGS_GATK_IND/$FILE_NAME.stderr

        echo "Processing of the file $FILE_NAME is finished :-) !!!"
        echo "##################################################"

done

echo "GATK (variant calling step (default mode)) ---- FINISHED"

# Fourth part: GATK based variant-calling [GVCF variant-calling]

echo "GATK (variant calling step (GVCF mode)) ---- STARTING"

for BAM in $OUT_SAMTOOLS/*.marked_duplicates.bam

    do

        FILE_NAME=$(basename $BAM | sed "s/.marked_duplicates.bam//" - )

        echo "##################################################"

        date
        echo "Sample: $FILE_NAME"

        echo "##################################################"

        time  gatk --java-options "-Xmx4g" HaplotypeCaller  \
                -R $RESOURCES/Cpepo_genome_v4.1.fa \
                -I $BAM \
                -O $OUT_GATK_JOINT/$FILE_NAME.g.vcf.gz \
                -ERC GVCF > $LOGS_GATK_JOINT/$FILE_NAME.stdout  2> $LOGS_GATK_JOINT/$FILE_NAME.stderr

        echo "Processing of the file $FILE_NAME is finished :-) !!!"
        echo "##################################################"

done

echo "GATK (variant calling step (GVCF mode)) ---- FINISHED"