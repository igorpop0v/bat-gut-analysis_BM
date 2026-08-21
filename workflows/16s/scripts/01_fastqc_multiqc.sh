#!/usr/bin/env bash

# Perform quality control of raw 16S rRNA FASTQ reads.
# FastQC analyzes each FASTQ file separately.
# MultiQC combines all FastQC reports into a single HTML report.
#
# Usage:
#   ./01_fastqc_multiqc.sh INPUT_DIR OUTPUT_DIR [THREADS]

set -euo pipefail

# Version-pinned Docker images ensure reproducibility.
readonly FASTQC_IMAGE="quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
readonly MULTIQC_IMAGE="multiqc/multiqc:v1.35"

# Determine the repository root from the script location.
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

if (( $# < 2 || $# > 3 )); then
    echo "Usage:"
    echo "  $0 INPUT_DIR OUTPUT_DIR [THREADS]"
    echo
    echo "Example:"
    echo "  $0 /path/to/raw_fastq /path/to/qc_results 8"
    exit 1
fi

RAW_DIR="$1"
QC_DIR="$2"
THREADS="${3:-8}"

# Validate command-line arguments.
if [[ ! -d "${RAW_DIR}" ]]; then
    echo "ERROR: input directory does not exist: ${RAW_DIR}" >&2
    exit 1
fi

if [[ ! "${THREADS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: THREADS must be a positive integer." >&2
    exit 1
fi

command -v docker >/dev/null 2>&1 || {
    echo "ERROR: Docker is not available." >&2
    exit 1
}

RAW_DIR="$(realpath "${RAW_DIR}")"
QC_DIR="$(realpath -m "${QC_DIR}")"

# Raw data and generated results must remain outside the Git repository.
if [[ "${RAW_DIR}" == "${REPO_ROOT}" ||
      "${RAW_DIR}" == "${REPO_ROOT}/"* ]]; then
    echo "ERROR: raw FASTQ files must be outside the Git repository." >&2
    exit 1
fi

if [[ "${QC_DIR}" == "${REPO_ROOT}" ||
      "${QC_DIR}" == "${REPO_ROOT}/"* ]]; then
    echo "ERROR: QC results must be written outside the Git repository." >&2
    exit 1
fi

readonly FASTQC_DIR="${QC_DIR}/fastqc"
readonly MULTIQC_DIR="${QC_DIR}/multiqc"
readonly LOG_DIR="${QC_DIR}/logs"

mkdir -p "${FASTQC_DIR}" "${MULTIQC_DIR}" "${LOG_DIR}"

# Find FASTQ files located directly in the input directory.
# Both forward (R1) and reverse (R2) reads must be included.
FASTQ_FILES=()

mapfile -d '' FASTQ_FILES < <(
    find "${RAW_DIR}" -maxdepth 1 -type f \
        \( -iname "*.fastq" \
           -o -iname "*.fastq.gz" \
           -o -iname "*.fq" \
           -o -iname "*.fq.gz" \) \
        -print0 |
        sort -z
)

if (( ${#FASTQ_FILES[@]} == 0 )); then
    echo "ERROR: no FASTQ files were found in: ${RAW_DIR}" >&2
    exit 1
fi

echo "FASTQ files found: ${#FASTQ_FILES[@]}"
echo "Input directory: ${RAW_DIR}"
echo "Output directory: ${QC_DIR}"
echo "FastQC threads: ${THREADS}"

# Convert host paths to paths available inside the Docker container.
CONTAINER_FASTQ_FILES=()

for fastq_file in "${FASTQ_FILES[@]}"; do
    filename="$(basename "${fastq_file}")"
    CONTAINER_FASTQ_FILES+=("/data/${filename}")
done

echo
echo "Running FastQC..."

docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --volume "${RAW_DIR}:/data:ro" \
    --volume "${FASTQC_DIR}:/output" \
    "${FASTQC_IMAGE}" \
    fastqc \
        --threads "${THREADS}" \
        --format fastq \
        --outdir /output \
        "${CONTAINER_FASTQ_FILES[@]}" \
    2>&1 | tee "${LOG_DIR}/fastqc.log"

echo
echo "Running MultiQC..."

docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --volume "${FASTQC_DIR}:/input:ro" \
    --volume "${MULTIQC_DIR}:/output" \
    "${MULTIQC_IMAGE}" \
    multiqc /input \
        --outdir /output \
        --filename multiqc_report.html \
        --force \
    2>&1 | tee "${LOG_DIR}/multiqc.log"

echo
echo "Quality control completed."
echo "MultiQC report:"
echo "${MULTIQC_DIR}/multiqc_report.html"
