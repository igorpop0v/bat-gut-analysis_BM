#!/usr/bin/env bash

set -e

VIROME_DIR="$1"

REFERENCE_DIR="${VIROME_DIR}/host_reference"
DOWNLOAD_DIR="${REFERENCE_DIR}/downloads"
INDEX_DIR="${REFERENCE_DIR}/bowtie2_index"

BOWTIE2_IMAGE="quay.io/biocontainers/bowtie2:2.5.5--ha27dd3b_0"

mkdir -p "$DOWNLOAD_DIR"
mkdir -p "$INDEX_DIR"

echo "Downloading Nyctalus leisleri genome..."

wget -c \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/264/875/GCA_964264875.2_mNycLei1.hap1.2/GCA_964264875.2_mNycLei1.hap1.2_genomic.fna.gz" \
  -O "${DOWNLOAD_DIR}/nyctalus_leisleri.fna.gz"

echo "Downloading Nyctalus aviator genome..."

wget -c \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/036/971/965/GCA_036971965.1_NENU_Navi_1.0/GCA_036971965.1_NENU_Navi_1.0_genomic.fna.gz" \
  -O "${DOWNLOAD_DIR}/nyctalus_aviator.fna.gz"

echo "Downloading Nyctalus noctula mitochondrial genome..."

wget \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_027237.1&rettype=fasta&retmode=text" \
  -O "${DOWNLOAD_DIR}/nyctalus_noctula_mitochondrion.fasta"

echo "Combining host reference sequences..."

gzip -dc "${DOWNLOAD_DIR}/nyctalus_leisleri.fna.gz" \
  > "${REFERENCE_DIR}/nyctalus_host_combined.fasta"

gzip -dc "${DOWNLOAD_DIR}/nyctalus_aviator.fna.gz" \
  >> "${REFERENCE_DIR}/nyctalus_host_combined.fasta"

cat "${DOWNLOAD_DIR}/nyctalus_noctula_mitochondrion.fasta" \
  >> "${REFERENCE_DIR}/nyctalus_host_combined.fasta"

echo "Building Bowtie2 index..."

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --volume "${REFERENCE_DIR}:/reference" \
  "$BOWTIE2_IMAGE" \
  bowtie2-build \
  --threads 12 \
  --large-index \
  /reference/nyctalus_host_combined.fasta \
  /reference/bowtie2_index/nyctalus_host

echo "Host reference preparation completed."
echo "Reference: ${REFERENCE_DIR}/nyctalus_host_combined.fasta"
echo "Bowtie2 index: ${INDEX_DIR}/nyctalus_host"
