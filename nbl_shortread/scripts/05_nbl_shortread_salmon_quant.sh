#!/bin/bash 
#SBATCH --job-name novel_tx_quant

# This Bash script performs parallel quantification of novel transcripts in RNA sequencing data using Salmon
# by parallel submission of a salmon quant bash script for each sample.

#########################################
#
# ADAPTED FOR NBL SHORT-READ DATA
#
#########################################

# Configuration (to be adjusted according to the dataset)
wd=`pwd`
sample_data_folder="${wd}/data/"
paired_end=true
index="${wd}/analysis/salmon_index/nbl_shortread_transcriptome"
transcriptome_basename="nbl_shortread_transcriptome"

# Default parameters
scriptdir="${wd}/scripts"
data_folder="/hpc/pmc_vanheesch/data"  # Data folder containing all of our sequencing data
simul_array_runs=50

run=$(echo 1-`date +%F`)

# Functions
function info() {
    echo "INFO: $@" >&2
}
function error() {
    echo "ERR:  $@" >&2
}
function fatal() {
    echo "ERR:  $@" >&2
    exit 1
}
function get_samples() {
  project_data_folder=$1
  data_folder=$2
  paired_end=${3:-true}
  
  # Check whether the files are in correct format
  if [[ $(find ${project_data_folder} -name "*.fastq.gz" | wc -l) -eq 0 ]]; then
    fatal "No .fastq.gz files found in ${project_data_folder} or its subfolders"
  else
    # Find unique R1 filenames and get corresponding R1/R2 fastq files
    r1_files=()
    r2_files=()
    if [[ "${paired_end,,}" == "true" ]]; then
      readarray -t r1_filenames_raw < <(find "${project_data_folder}/" -maxdepth 3 -regextype posix-extended -regex '.*_R1(_001)?\.fastq\.gz$' -printf '%f\n' | sort -u)
      for r1_filename in "${r1_filenames_raw[@]}"; do
        full_path="${project_data_folder}/${r1_filename}"
        # If r1_filename is a symlink, find original file
        if [[ -L "${full_path}" ]]; then
          resolved_path="$(readlink -f "${full_path}")"
          if [[ ! -e "${resolved_path}" ]]; then
            echo "Warning: Broken symlink for ${r1_filename}. Skipping..."
            continue
          else
            r1_file="${resolved_path}"
          fi
        elif [[ -f "${full_path}" ]]; then
          r1_file="${full_path}"
        else
          echo "Warning: File ${r1_filename} not found. Skipping..."
          continue
        fi

        r2_file="$(echo "${r1_file}" | sed 's/_R1\(_001\)\?\.fastq\.gz/_R2\1.fastq.gz/')"
        if [[ ! -f "${r2_file}" ]]; then

          fatal "R2 file ${r2_file} not found for ${r1_file}"
        fi
        
        r1_files+=("${r1_file}")
        r2_files+=("${r2_file}")
      done
    else
      readarray -t r1_filenames < <(find "${project_data_folder}/" -maxdepth 1 -name "*_R1*" -printf '%f\n' | sort -u)
      for r1_filename in "${r1_filenames[@]}"; do
        if [[ -L "${project_data_folder}/${r1_filename}" ]]; then
            if [[ -e "${project_data_folder}/${r1_filename}" ]]; then
              r1_file="$(readlink -f "${project_data_folder}/${r1_filename}")"
            else
              echo "Warning: Broken symlink for ${r1_filename}. Skipping..."
              continue
            fi
        elif [[ -f "${project_data_folder}/${r1_filename}" ]]; then
            r1_file="${project_data_folder}/${r1_filename}"
        else
            echo "Warning: File ${r1_filename} not found. Skipping..."
            continue
        fi

        r1_files+=("${r1_file}")
      done
    fi
  fi

  # Initiate arrays
  sample_ids=()
  samples=()

  # Get sample IDs from fastq files
  for r1_file in "${r1_files[@]}"; do
    sample=$(basename "${r1_file}")
    #sample_id=$(basename ${r1_file} | rev | cut -d '_' -f 2- | rev | sort | uniq)
    sample_id=$(basename ${r1_file} | sed 's/_R1\(_[0-9]*\)\?.fastq.gz$//')
    samples+=("${sample}")
    sample_ids+=("${sample_id}")
  done

  # Make sure there are samples
  if [[ ${#samples[@]} -eq 0 ]]; then
    fatal "No samples found in ${project_data_folder}/"
  fi

  info "Samples:"
  for i in ${!samples[@]}; do
    info "$((i+1))    ${samples[i]}"
  done

  export r1_files
  export r2_files
  export sample_ids
  export samples
}

# get_samples $sample_data_folder $data_folder $paired_end

# printf "%s\n" "${r1_files[@]}" > ${wd}/documentation/r1_files.txt
# printf "%s\n" "${r2_files[@]}" > ${wd}/documentation/r2_files.txt
# printf "%s\n" "${sample_ids[@]}" > ${wd}/documentation/sample_ids.txt

echo "`date` running salmon quant run $run"

mkdir -p "log/${run}/${transcriptome_basename}/"

mapfile -t sample_ids < "${wd}/documentation/sample_ids.txt"

salmon_quant_jobid=()

salmon_quant_jobid+=($(sbatch --parsable \
  --mem=48G \
  --cpus-per-task=6 \
  --time=36:00:00 \
  --gres=tmpspace:50G \
  --array 1-${#sample_ids[@]}%${simul_array_runs} \
  --job-name=${run}.quant_novel_txome \
  --output=log/${run}/%A_%a.out \
  --export=ALL \
  ${scriptdir}/run_novel_tx_quant.sh
))

echo "Novel salmon jobid: ${salmon_quant_jobid[@]}"

#,wd=${wd},index=${index},transcriptome_basename=${transcriptome_basename} \