#!/bin/bash -e
#SBATCH -p jic-short
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=judit.talas@jic.ac.uk
#SBATCH --array=0-3
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=bsseq-trim

source cutadapt-1.9.1
source trim_galore-0.4.2

libraries=(JWXF001C_SDXF021 JWXF001C_SDXF022 JWXF001D_SDXF021 JWXF001D_SDXF022)

trim_galore ${libraries[$SLURM_ARRAY_TASK_ID]}.fastq.gz



































#srun bismark_methylation_extractor --comprehensive --cytosine_report --bedGraph --CX ${ARRAY[$SLURM_ARRAY_TASK_ID]}.bam --genome_folder ./ref2 --single-end --gazillion
#cat ${ARRAY[$SLURM_ARRAY_TASK_ID]}.CX_report.txt | awk '{if($6=="CHH") print $0}' | awk '{if($7=="CAT" || $7=="CAA" || $7=="CAC") print $0}' > ${ARRAY[$SLURM_ARRAY_TASK_ID]}.CAH_CX_report.txt
#cat ${ARRAY[$SLURM_ARRAY_TASK_ID]}.CX_report.txt | awk '{if($6=="CHH") print $0}' | awk '{if($7=="CCT" || $7=="CCA" || $7=="CCC") print $0}' > ${ARRAY[$SLURM_ARRAY_TASK_ID]}.CCH_CX_report.txt
#cat ${ARRAY[$SLURM_ARRAY_TASK_ID]}.CX_report.txt | awk '{if($6=="CHH") print $0}' | awk '{if($7=="CTT" || $7=="CTA" || $7=="CTC") print $0}' > ${ARRAY[$SLURM_ARRAY_TASK_ID]}.CTH_CX_report.txt

#perl -S ~/groupdir/bin/cytosine_report2gff.pl ${ARRAY[$SLURM_ARRAY_TASK_ID]}.CX_report.txt
#perl -S ~/groupdir/bin/cytosine_report2gff.pl ${ARRAY[$SLURM_ARRAY_TASK_ID]}.CAH_CX_report.txt
#perl -S ~/groupdir/bin/cytosine_report2gff.pl ${ARRAY[$SLURM_ARRAY_TASK_ID]}.CCH_CX_report.txt
#perl -S ~/groupdir/bin/cytosine_report2gff.pl ${ARRAY[$SLURM_ARRAY_TASK_ID]}.CTH_CX_report.txt


#~/groupdir/Martin_Vickers/bin/nightfury -i ${ARRAY[$SLURM_ARRAY_TASK_ID]}.fq.gz -r MpTak1v5.1_and_plastids.fasta -nd -o ${ARRAY[$SLURM_ARRAY_TASK_ID]}.nd.bam
#srun bismark_methylation_extractor --comprehensive --cytosine_report --bedGraph --CX ${ARRAY[$SLURM_ARRAY_TASK_ID]}.nd.bam --genome_folder ./ref2 --single--end --gazillion
#perl -S ~/groupdir/bin/cytosine_report2gff.pl ${ARRAY[$SLURM_ARRAY_TASK_ID]}.nd.CX_report.txt

