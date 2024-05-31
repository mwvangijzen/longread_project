#!/bin/bash
#SBATCH --job-name=nfMerge
#SBATCH --output=/hpc/pmc_vanheesch/projects/mvangijzen/nbl_longread/log/nfMerge%A.out
#SBATCH --error=/hpc/pmc_vanheesch/projects/mvangijzen/nbl_longread/log/nfMerge%A.err
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
#SBATCH --nodes=1
#SBATCH --open-mode=append



module load java
module load nextflow/23.04.4

nextflow run mars13/nf_rna_pipeline \
    -c documentation/gffcompare_params.config \
    -profile slurm