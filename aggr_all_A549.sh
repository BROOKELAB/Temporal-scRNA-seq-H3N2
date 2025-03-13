#!/bin/bash
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=150G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=joel
#SBATCH -J aggr-all
#SBATCH --output=aggr-all-%j.out
#SBATCH -D ../slurm_all


### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

### Load Modules
module load cellranger/7.1.0

### specifying the output folder
cd ../results/cellranger

cellranger aggr --id=Aggregated_samples \
                  --csv=../../src/AggrList081524.csv \
                  --normalize=none
