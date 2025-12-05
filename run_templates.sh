#!/bin/bash
#SBATCH -A b1042             ## account (b1012 for buyin, genomics = b1042)
#SBATCH -p genomics          ## "-p" instead of "-q"
#SBATCH -J run_seurat_pipeline        ## job name
#SBATCH --mail-type=FAIL,TIME_LIMIT_90
#SBATCH --mail-user=your_email@northwestern.edu
#SBATCH -o "%x.o%j"
#SBATCH -N 1                 ## number of nodes
#SBATCH -n 1                 ## number of cores
#SBATCH -t 48:00:00          ## walltime
#SBATCH --mem=80G

export MC_CORES=${SLURM_NTASKS}

module purge all
module load R/4.4.0
module load geos/3.8.1
module load hdf5/1.8.19-serial
module load pandoc/2.2.1

number=$1
batch=$2

echo -n  "Running pipeline with arguments report number = ${number}, batch control method = ${batch}, at "
date

SECONDS=0

# Rscript ${script}
Rscript execute_pipeline-Ex.R ${number} ${batch}

duration=$SECONDS
printf "\nProcessed in %.2f minutes\n" $(($duration / 60))
