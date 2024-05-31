#!/bin/bash 
wd=`pwd`
index="${wd}/analysis/salmon_index/nbl_shortread_transcriptome"
transcriptome_basename="nbl_shortread_transcriptome"
data_folder="/hpc/pmc_vanheesch/data"
cpu=6

# Load correct modules
module load salmon/1.8.0
# module load cutadapt/3.4
# module load fastqc/0.11.9
# module load trimgalore/0.6.6

# Get correct files
mapfile -t r1_files < "${wd}/documentation/r1_files.txt"
mapfile -t r2_files < "${wd}/documentation/r2_files.txt"
mapfile -t sample_ids < "${wd}/documentation/sample_ids.txt"

full_path_fastq_1="${r1_files[$((SLURM_ARRAY_TASK_ID-1))]}"
full_path_fastq_2="${r2_files[$((SLURM_ARRAY_TASK_ID-1))]}"

# Set names
sample_id="${sample_ids[$((SLURM_ARRAY_TASK_ID-1))]}"

bf1=$(basename ${full_path_fastq_1})
bf2=$(basename ${full_path_fastq_2})

r1_trimmed="${bf1/_R1_/_R1_trimmed_}"
r2_trimmed="${bf2/_R2_/_R2_trimmed_}"

# Check if "quant.sf" file in salmon output dir already exists for this sample, if so, skip this sample
if [[ -f "${wd}/analysis/salmon_quant/${sample_id}/quant.sf" ]]; then
  echo "quant.sf file already exists for ${sample_id}, skipping..."
  exit 0
fi

# Create output dirs
# echo "`date` running trimgalore for ${sample_id}"

# Trim reads with fastp
apptainer exec -B "/hpc:/hpc","${TMPDIR}:${TMPDIR}" --env LC_ALL=C.UTF-8 /hpc/local/Rocky8/pmc_vanheesch/singularity_images/fastp_0.23.4--hadf994f_2.sif fastp \
  -i "${full_path_fastq_1}" \
  -I "${full_path_fastq_2}" \
  -o "${TMPDIR}/${bf1}" \
  -O "${TMPDIR}/${bf2}" \
  --detect_adapter_for_pe \
  --verbose \
  --thread "${cpu}"

# # Run trimgalore on both reads

# cd "${TMPDIR}"

# trim_galore "${full_path_fastq_1}" "${full_path_fastq_2}" \
#   --cores 2 \
#   --paired \
#   --gzip \
#   --output_dir "${TMPDIR}"

# # Change names of validated trimgalore output to basename
# mv "${TMPDIR}/${bf1%.*.*}_val_1.fq.gz" "${TMPDIR}/${bf1}"
# mv "${TMPDIR}/${bf2%.*.*}_val_2.fq.gz" "${TMPDIR}/${bf2}"

# echo "`date` Running salmon quant for ${sample_id}"

echo "`date` Running salmon quant with ${transcriptome_basename} for ${sample_id}"

mkdir -p "${wd}/analysis/salmon_quant/${sample_id}/"
# Run salmon for transcript counts

#apptainer exec -B "/hpc:/hpc","${TMPDIR}:${TMPDIR}"  --env LC_ALL=C.UTF-8 /hpc/local/Rocky8/pmc_vanheesch/singularity_images/salmon-1.8.0.sif salmon quant \

salmon quant \
  --libType "A" \
  --validateMappings \
  --gcBias \
  --numGibbsSamples 30 \
  --threads "${cpu}" \
  -i "${index}" \
  -1 "${TMPDIR}/${bf1}" \
  -2 "${TMPDIR}/${bf2}" \
  --output "${wd}/analysis/salmon_quant/${sample_id}/"

echo "`date` Finished salmon ${sample_id}"
