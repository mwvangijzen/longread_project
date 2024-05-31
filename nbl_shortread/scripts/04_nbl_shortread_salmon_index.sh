#!/bin/bash
#SBATCH --job-name=salmon_index
#SBATCH --output=log/salmon_index%A.out
#SBATCH --error=log/salmon_index%A.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --export=NONE
#SBATCH --gres=tmpspace:20G

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date:    24-06-2022
# 
# Adapted: MG 
# Date:    28-02-2024
#
######################################################################

set -uo pipefail

# Load parameters from main script
cpu=8
gffread_version=0.12.7
salmon_version=1.8.0
wd=`pwd`
resource_dir="/hpc/pmc_vanheesch/shared_resources"
species="Homo_sapiens"
genome_version="GRCh38"
annot="102"
scriptdir="${wd}/scripts"
reference_genome="/${resource_dir}/GENOMES/${species}.${genome_version}/${annot}/${species}.${genome_version}.dna.primary_assembly.fa"

merged_gtf_basename="nbl_shortread_transcriptome"

#from orignal script: merged_gtf_basename=${gtf_basenames[$((SLURM_ARRAY_TASK_ID-1))]}

# Load correct modules
# module load gffread/${gffread_version}
module load salmon/${salmon_version}

# Create output dirs
mkdir -p "${wd}/analysis/salmon_index/"

# Create fasta from merged annotated GTF
apptainer exec -B "/hpc:/hpc" --env LC_ALL=C.UTF-8 /hpc/local/Rocky8/pmc_vanheesch/singularity_images/gffread-0.12.6.sif gffread \
  -w "${wd}/analysis/salmon_index/${merged_gtf_basename}_transcripts.fa" \
  -g "${reference_genome}" \
  "${wd}/analysis/customannotation/${merged_gtf_basename}_novel_filtered_new_new_len3.gtf"

# Create "gentrome" for salmon index
gentrome="${wd}/analysis/salmon_index/gentrome_${merged_gtf_basename}.fa"
cat "${wd}/analysis/salmon_index/${merged_gtf_basename}_transcripts.fa" "${reference_genome}" > "${gentrome}"

# Grep seq names from fa and create decoys for salmon
decoy="${wd}/analysis/salmon_index/decoys_${merged_gtf_basename}.txt"
grep "^>" "${reference_genome}" | cut -d " " -f 1 > "${decoy}"
sed -i.bak -e 's/>//g' "${decoy}"

# Create salmon index
salmon index \
  --transcripts "${gentrome}" \
  --decoys "${decoy}" \
  --index "${wd}/analysis/salmon_index/${merged_gtf_basename}" \
  --threads "${cpu}"

echo "`date` Finished creating salmon index"