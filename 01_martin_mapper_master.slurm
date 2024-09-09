#!/bin/bash -e
#SBATCH -p jic-medium
#SBATCH --mem 40G
#SBATCH -c 16
#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=judit.talas@jic.ac.uk
#SBATCH --job-name=bs_map
#SBATCH --array=0

libraries=(Mp_WT_embryo_merged Mp_WT_thallus_merged)

work_folder="/jic/scratch/groups/Phil-Carella/talas/Marchantia/BS-seq"
genomes="/jic/scratch/groups/Phil-Carella/talas/Marchantia/BS-seq/annotations"

##run the below for the first time on just one library as it creates mapped intermediate files. Legacy methods are included in original script found in the following directory:
##/jic/scratch/groups/Xiaoqi-Feng/talas/Scripts/Jincheng_Jimmy/BS_Seq/

# map to genome

~/groupdir/Martin_Vickers/bin/nightfury -i $work_folder/${libraries[$SLURM_ARRAY_TASK_ID]}_trimmed.fq.gz -r $genomes/tak_v6/MpTak_v6.1.genome.autosomes.fasta -o $work_folder/takv6_autosomes/${libraries[$SLURM_ARRAY_TASK_ID]}.bam

# sort and index bam files

cd $work_folder/takv6_autosomes

#sometimes samtools fails in this script - in that case run samtools commands interactively then continue

source samtools-1.9

samtools sort -o sorted_${libraries[$SLURM_ARRAY_TASK_ID]}.bam ${libraries[$SLURM_ARRAY_TASK_ID]}.bam 
samtools index sorted_${libraries[$SLURM_ARRAY_TASK_ID]}.bam


source MethylDackel-0.4.0
#map using MethylDackel and extract non CG reads - results in 3 bedGraph files 
srun MethylDackel extract \
--CHG \
--CHH \
-@ 15 \
$genomes/tak_v6/MpTak_v6.1.genome.autosomes.fasta \
sorted_${libraries[$SLURM_ARRAY_TASK_ID]}.bam


# collapse to w100 and edit label field to say CG CHG and CHH for SignalMap compatibility

for methylation_context in CpG CHG CHH
    do
        cat ${libraries[$SLURM_ARRAY_TASK_ID]}_"$methylation_context".bedGraph | \
            grep -v "track" | \
            awk '{print $1 "\t.\t.\t"$2+1 "\t"$3+1 "\t"$4/100"\t.\t.\tc="$5";t="$6}' \
             >  ${libraries[$SLURM_ARRAY_TASK_ID]}_"$methylation_context".w1.gff
        if [ "$methylation_context" == "CpG" ] ; then
                ~/groupdir/Martin_Vickers/bin/w50_creator \
                -i  ${libraries[$SLURM_ARRAY_TASK_ID]}_"$methylation_context".w1.gff \
                -l  ${libraries[$SLURM_ARRAY_TASK_ID]}_CG -s 100 \
                > ${libraries[$SLURM_ARRAY_TASK_ID]}."$methylation_context".w100.gff
        else
                ~/groupdir/Martin_Vickers/bin/w50_creator \
                -i  ${libraries[$SLURM_ARRAY_TASK_ID]}_"$methylation_context".w1.gff \
                -l  ${libraries[$SLURM_ARRAY_TASK_ID]}_"$methylation_context" -s 100 \
               > ${libraries[$SLURM_ARRAY_TASK_ID]}."$methylation_context".w100.gff
        fi
    done