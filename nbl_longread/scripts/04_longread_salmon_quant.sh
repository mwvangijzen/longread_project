#!/bin/bash
#SBATCH --job-name=longread_salmon_quant
#SBATCH --output=log/lr_salmon_quant%A.out
#SBATCH --error=log/lr_salmon_quant%A.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --export=NONE
#SBATCH --gres=tmpspace:250G
#SBATCH --array 1-6%6

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

# Set names
SAMPLE_IDS=("A13" "A14" "A15" "E02" "E03" "E05")
S="${SAMPLE_IDS[$((SLURM_ARRAY_TASK_ID-1))]}"
FASTQ_IDS=("A13-047-cells-TeloP-TotalRNA_1.fastq.gz" "A14-048-cells-TeloP-TotalRNA_1.fastq.gz" "A15-049-cells-TeloP-TotalRNA_1.fastq.gz" \
"E02-036-cells-TeloP-TotalRNA_1.fastq.gz" "E03-037-cells-TeloP-TotalRNA_1.fastq.gz" "E05-039-cells-TeloP-TotalRNA_1.fastq.gz")
FASTQ="${FASTQ_IDS[$((SLURM_ARRAY_TASK_ID-1))]}"

# Align reads to the transcriptome
## Perform alignment
apptainer exec -B /hpc:/hpc,${TMPDIR}:${TMPDIR} "${MINIMAP_CONTAINER}" minimap2 -ax map-ont \
"${WD}/analysis/salmon_count/nbl_longread_transcriptome_transcripts.fa" \
"${WD}/data/${S}/${FASTQ}" > \
"${TMPDIR}/${S}_transcriptome_aln.sam" 

## Convert sam to bam
apptainer exec -B /hpc:/hpc,${TMPDIR}:${TMPDIR} "${SAMTOOLS_CONTAINER}" samtools view -S -b \
"${TMPDIR}/${S}_transcriptome_aln.sam" > \
"${WD}/analysis/salmon_count/transcriptome_alignment/${S}_transcriptome_aln.bam"

# Alignment-based quantification
apptainer exec -B /hpc:/hpc,${TMPDIR}:${TMPDIR} "${SALMON_CONTAINER}" salmon quant \
    --noErrorModel \
    -p "${CPU}" \
    -t "${REF_TRANSCRIPTOME}" \
    -l SF \
    -a "${BAM_DIR}/${S}_transcriptome_aln.bam" \
    -o "${TMPDIR}"
    mv "${TMPDIR}/quant.sf" "${WD}/analysis/salmon_count/counts/${S}_${MERGED_GTF_BASENAME}.transcript_counts.tsv"

# Can you only make it echo this sentence after it finishes the other steps?
echo "`date` Finished creating salmon counts"