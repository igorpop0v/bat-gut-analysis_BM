# 16S rRNA read-quality assessment and processing decision

## Dataset

- Target: V3–V4 region of the bacterial 16S rRNA gene
- Sequencing format: paired-end, 2 × 301 bp
- Number of samples: 77
- Library preparation kit: Quick-16S NGS Library Prep Kit
  (Zymo Research, Irvine, CA, USA)

## Primers

Two variants of the 341F primer were used:

- `CCTACGGGDGGCWGCAG`
- `CCTAYGGGGYGCWGCAG`

Reverse primer 806R:

- `GACTACNVGGGTMTCTAATCC`

## Quality-control software

- FastQC 0.12.1
- MultiQC 1.35

Raw reads and complete quality-control reports are stored outside this
repository.

## Quality-control observations

FastQC reports were generated separately for all forward and reverse read
files and summarized with MultiQC.

The forward reads retained high per-base quality over a substantially longer
part of the read than the reverse reads.

Approximate mean Phred quality scores across samples:

| Position | Forward reads (R1) | Reverse reads (R2) |
|---:|---:|---:|
| 150 | 36.5 | 32.0 |
| 180 | 35.9 | 28.7 |
| 200 | 35.7 | 25.0 |
| 220 | 34.4 | 18.5 |
| 230 | 32.7 | 14.1 |
| 250 | 29.8 | 11.2 |
| 300 | 16.3 | 10.2 |

Reverse-read quality generally decreased below Q30 at approximately
170 nucleotides and below Q20 at approximately 220 nucleotides.

Forward-read quality generally remained above Q30 until approximately
230–240 nucleotides.

No substantial adapter-content problem was detected.

High sequence duplication, non-random base composition, and overrepresented
sequences were not used as exclusion criteria because these patterns are
expected for targeted 16S rRNA amplicon sequencing.

## Paired-end overlap assessment

The expected length of the V3–V4 amplicon is approximately 460 nucleotides,
or approximately 420–430 nucleotides after primer removal.

After removing the low-quality portions of the reverse reads, the retained
forward and reverse sequences are not expected to provide sufficient
high-quality overlap for reliable paired-end merging.

Retaining additional reverse-read positions solely to obtain overlap would
include bases with low Phred quality and would be expected to increase read
loss during filtering and merging.

## Processing decision

The primary analysis will therefore use forward reads only.

The following preprocessing strategy will be applied:

1. Import forward reads as single-end Phred33 data.
2. Remove both variants of the 341F primer with Cutadapt.
3. Discard reads in which a forward primer cannot be identified.
4. Summarize primer-trimmed reads with QIIME 2.
5. Apply QIIME 2 quality filtering.
6. Denoise the filtered reads with Deblur.
7. Truncate primer-free reads to 190 nucleotides during Deblur denoising.
8. Generate per-sample denoising statistics.

The Deblur trim length of 190 nucleotides was selected conservatively. At the
corresponding position in the raw forward reads, the average Phred quality
remains high.

## Reproducibility

Software versions, primer sequences, processing parameters, and analysis
commands are maintained in this repository. Raw sequencing reads, QIIME 2
artifacts, and generated reports are excluded from version control.
