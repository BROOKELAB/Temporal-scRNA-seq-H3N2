#!/bin/bash
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=150G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=joel
#SBATCH -J CR-count-perth09-timepoint
#SBATCH --output=cell-ranger-%A_%a.out
#SBATCH --array=1-11%1
#SBATCH -D ../slurm_all


cd ../

Sample_ID=`head files_perth09.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 1`
Sample_name=`head files_perth09.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 2`
File_folder=`head files_perth09.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 3`

### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

### Load Modules
module load cellranger/7.1.0

### specifying the output folder

mkdir -p ../results/cellranger
cd ../results/cellranger


# Try with ref made with minimal GTF:
DB=../cellranger_4.0.0_hg38_plus_perth1609_renamed



### Run app on file
cellranger count --id=$Sample_ID \
    --localcores=$SLURM_CPUS_ON_NODE  \
    --localmem 128 \
    --transcriptome=$DB \
    --fastqs=$File_folder \
    --sample=$Sample_name \
    --expect-cells=5000

