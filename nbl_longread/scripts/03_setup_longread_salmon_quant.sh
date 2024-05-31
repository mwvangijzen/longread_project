#!/bin/bash
#SBATCH --job-name=longread_salmon_setup
#SBATCH --output=log/lr_salmon_setup%A.out
#SBATCH --error=log/lr_salmon_setup%A.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --export=NONE
#SBATCH --gres=tmpspace:20G

######################################################################
# Author:  MG
# Date:    29-02-2024
######################################################################

set -uo pipefail

# Parse command-line arguments
while [[ $# -gt 0 ]]
do
key="$1"
 
case $key in
    -c|--config)
    CONFIG="$2"
    shift
    shift
    ;;
    -h|--help)
    echo "See Marina"
    exit 1
    ;;
    "")
    echo "Error: no option provided"
    exit 1
    ;;
    *)
    echo "Unknown option: $key"
    exit 1
    ;;
esac
done

#Alternatively provide a config file with all paths and params with:
source "${CONFIG}"

# Create output dirs
mkdir -p "${WD}/analysis/salmon_count/"

# Create fasta from merged annotated GTF
apptainer exec -B /hpc:/hpc "${GFFREAD_CONTAINER}" gffread \
  -w "${WD}/analysis/salmon_count/${MERGED_GTF_BASENAME}_transcripts.fa" \
  -g "${REF_GENOME}" \
  "${WD}/analysis/gffcompare_combine/customannotation/${MERGED_GTF_BASENAME}_novel_filtered_fixed.gtf"

# Create minimap2 index
apptainer exec -B /hpc:/hpc "${MINIMAP_CONTAINER}" minimap2 \
 -d "${WD}/analysis/salmon_count/transcriptome_alignment/longread_transcriptome.mmi" \
 "${WD}/analysis/salmon_count/nbl_longread_transcriptome_transcripts.fa"