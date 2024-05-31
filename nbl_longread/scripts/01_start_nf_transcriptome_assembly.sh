#!/bin/bash
#SBATCH --job-name=nfTXome
#SBATCH --output=log/nfTXome%A.out
#SBATCH --error=log/nfTXome%A.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --export=NONE


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

#Loads default java version (if not 17, specify)
module load java
#Loads default nextflow version (if not >=23.04.2, specify)
module load nextflow/23.04.4

#Params with -- belong to the workflow, params with - belong to nextflow
#nextflow profile defaults to docker, need to specify singularity

nextflow run epi2me-labs/wf-transcriptomes \
  --fastq ${FASTQ} \
  --ref_genome ${REF_GENOME} \
  --ref_annotation ${REF_ANNOT} \
  --out_dir ${OUTDIR} \
  --sample_sheet ${SAMPLE_SHEET} \
  --minimap2_index_opts '-k 15' \
  --stringtie_opts '-s 1000' \
  --bundle_min_reads '50000' \
  -profile ${NF_PROFILE} \
  -c ${NF_CONFIG} \
  -with-report ${WD}/log/resources.html \
  -resume #use this option to resume the workflow after the last fully completed step from the previous run