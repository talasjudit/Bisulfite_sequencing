#!/bin/bash -e
#SBATCH -p nbi-medium
#SBATCH --mem=12000
#SBATCH --cpus-per-task=4
#SBATCH --job-name=mashup
#SBATCH --array=0-5

source bssequel-0.0.1
source xftools-0.0.1
source samtools-1.7

ARRAY=(
Mp_WT_embryo_merged_CpG \
Mp_WT_embryo_merged_CHG \
Mp_WT_embryo_merged_CHH \
Mp_WT_thallus_merged_CpG \
Mp_WT_thallus_merged_CHG \
Mp_WT_thallus_merged_CHH
)

# 5 primed and 3 primed ends analyisis for w1 and average ends
for primed in 3 5
    do
        srun /jic/research-groups/Xiaoqi-Feng/bin/ends_analysis.pl -g annotations/tak1_20190222_repeatmodeler_withdata.JW3.sorted.gff \
        -b 100 \
        -d 5000 \
        -s 6 \
        -x ID \
        -k 1500 \
        -"$primed" \
        ${ARRAY[$SLURM_ARRAY_TASK_ID]}.w1.gff > \
        ${ARRAY[$SLURM_ARRAY_TASK_ID]}.tak1_TE."$primed"ends

        # average ends
        srun /jic/research-groups/Xiaoqi-Feng/bin/average_ends_new.pl \
        -s 100 \
        -w 100 \
        ${ARRAY[$SLURM_ARRAY_TASK_ID]}.tak1_TE."$primed"ends > \ 
        ${ARRAY[$SLURM_ARRAY_TASK_ID]}.tak1_TE."$primed"ends.avg
    done