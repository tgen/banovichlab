#!/bin/bash
#SBATCH --job-name=velocyto
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=agutierrez@tgen.org
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 64000
#SBATCH --time=48:00:00
######################################################################
export PATH=$PATH:~/.local/bin

module load samtools/1.9

cd /scratch/agutierrez/10x_fastq/Outs/IPF/

velocyto run10x ${SAMPLE_DIR} ${REF_GTF}
