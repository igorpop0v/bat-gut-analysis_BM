# Software and container versions

The command-line workflows use fixed container tags wherever a stable versioned image was available.

| Analysis component | Software or image |
|---|---|
| Amplicon processing | `quay.io/qiime2/qiime2:2026.7` |
| Read QC | `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0` |
| QC aggregation | `multiqc/multiqc:v1.35` |
| Read trimming | `quay.io/biocontainers/fastp:1.3.6--h43da1c4_0` |
| Host-read removal | `quay.io/biocontainers/bowtie2:2.5.5--ha27dd3b_0` |
| rRNA removal | local image `bat-gut-sortmerna:7.0.0` built from `workflows/virome/containers/sortmerna/Dockerfile` |
| Read and contig classification | `ghcr.io/steineggerlab/metabuli:1.2.0` |
| Virome assembly | `quay.io/biocontainers/megahit:1.2.9--h5ca1c30_6` |
| Assembly QC | `staphb/quast:5.3.0` |
| Nucleotide and protein similarity searches | `ncbi/blast:2.17.0` |
| Densovirus read mapping | local image `bat-gut-densovirus-mapping:1.0` built from `workflows/virome/containers/densovirus_mapping/Dockerfile` |
| GenBank parsing | Python 3.12.11 and Biopython 1.85, built from `workflows/virome/containers/biopython/Dockerfile` |
| Reference protein clustering | `quay.io/biocontainers/proteinortho:6.3.6--h2b77389_0` |
| Protein alignment | `evolbioinfo/mafft:v7.526` |
| Alignment trimming | `evolbioinfo/trimal:v1.5.1` |
| Phylogenetic inference | `evolbioinfo/iqtree:v3.1.1` |
| Functional prediction | `quay.io/biocontainers/picrust2:2.6.3--pyhdfd78af_2` |

## Locally built images

Build the SortMeRNA image from the repository root:

```bash
docker build \
  --tag bat-gut-sortmerna:7.0.0 \
  workflows/virome/containers/sortmerna
```

Build the densovirus read-mapping image:

```bash
docker build \
  --tag bat-gut-densovirus-mapping:1.0 \
  workflows/virome/containers/densovirus_mapping
```

Build the Biopython image:

```bash
docker build \
  --tag bat-gut-biopython:1.85 \
  workflows/virome/containers/biopython
```

## Reference downloads

`workflows/virome/scripts/22_download_densovirus_reference_panel.sh` uses the NCBI EDirect container to retrieve the accession list recorded in `workflows/virome/config/densovirus_refseq_accessions_2024.txt`. The accession list, rather than a mutable search result, defines the reference panel used in the analysis.

The MetaBuli database is too large for Git. Its exact local snapshot and download date should be reported in the manuscript methods and retained with the project records.
