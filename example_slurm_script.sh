#!/bin/bash

#SBATCH --workdir=/work/m/maxwell9/
#SBATCH --job-name=seq07719
#SBATCH --time=01:55:00
#SBATCH --array=1-8
#SBATCH -N 1
#SBATCH -n 5
#SBATCH --dependency=afterok:3554649
#SBATCH --partition=development



reads[0]="/shares/biocomputing_i/FASTQ_Files/test2/PB58_36hrs_BP2_S3_L001_R1_001.fastq.gz"
reads[1]="/shares/biocomputing_i/FASTQ_Files/test2/PB58_36hrs_BP2_S3_L001_R2_001.fastq.gz"
reads[2]="/shares/biocomputing_i/FASTQ_Files/test2/PB58_12hrs_BP2_S1_L001_R1_001.fastq.gz"
reads[3]="/shares/biocomputing_i/FASTQ_Files/test2/PB58_12hrs_BP2_S1_L001_R2_001.fastq.gz"
reads[4]="/shares/biocomputing_i/FASTQ_Files/test2/PB58_24hrs_BP1_S6_L001_R1_001.fastq.gz"
reads[5]="/shares/biocomputing_i/FASTQ_Files/test2/PB58_24hrs_BP1_S6_L001_R2_001.fastq.gz"
reads[6]="/shares/biocomputing_i/FASTQ_Files/test2/PB58_48hrs_BP2_S8_L001_R1_001.fastq.gz"
reads[7]="/shares/biocomputing_i/FASTQ_Files/test2/PB58_48hrs_BP2_S8_L001_R2_001.fastq.gz"
reads[8]="/shares/biocomputing_i/FASTQ_Files/test2/PB58_36hrs_BP1_S7_L001_R1_001.fastq.gz"
reads[9]="/shares/biocomputing_i/FASTQ_Files/test2/PB58_36hrs_BP1_S7_L001_R2_001.fastq.gz"
reads[10]="/shares/biocomputing_i/FASTQ_Files/test2/PB58_48hrs_BP1_S4_L001_R1_001.fastq.gz"
reads[11]="/shares/biocomputing_i/FASTQ_Files/test2/PB58_48hrs_BP1_S4_L001_R2_001.fastq.gz"
reads[12]="/shares/biocomputing_i/FASTQ_Files/test2/PB58_24hrs_BP2_S2_L001_R1_001.fastq.gz"
reads[13]="/shares/biocomputing_i/FASTQ_Files/test2/PB58_24hrs_BP2_S2_L001_R2_001.fastq.gz"
reads[14]="/shares/biocomputing_i/FASTQ_Files/test2/PB58_12hrs_BP1_S5_L001_R1_001.fastq.gz"
reads[15]="/shares/biocomputing_i/FASTQ_Files/test2/PB58_12hrs_BP1_S5_L001_R2_001.fastq.gz"

echo "slurm job id: $SLURM_JOB_ID
slurm cluster name: $SLURM_CLUSTER_NAME
slurm nodelist: $SLURM_JOB_NODELIST"


tophat2 -G /shares/biocomputing_i/FASTQ_Files/test2/genes.gff -o thout_07719'_'$SLURM_ARRAY_TASK_ID -p 5 /work/m/maxwell9/genome07719 ${reads[$(($SLURM_ARRAY_TASK_ID * 2 - 2))]} ${reads[$(($SLURM_ARRAY_TASK_ID * 2 - 1))]}


cufflinks -p 5 -o clout_07719'_'$SLURM_ARRAY_TASK_ID thout_07719'_'$SLURM_ARRAY_TASK_ID/accepted_hits.bam


cuffquant -p 5 -o cqout_07719'_'$SLURM_ARRAY_TASK_ID /shares/biocomputing_i/FASTQ_Files/test2/genes.gff thout_07719'_'$SLURM_ARRAY_TASK_ID/accepted_hits.bam
