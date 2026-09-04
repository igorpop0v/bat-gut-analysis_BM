# Densovirus genome analysis

## Scope

This workflow evaluates densovirus-like sequences recovered from DNase-treated faecal RNA-sequencing libraries from common noctules (*Nyctalus noctula*). It combines read support, genome similarity, protein-marker phylogeny, and conservative feature annotation.

The presence of a sequence in a faecal RNA library is not treated as evidence of productive infection in bats. In insectivorous hosts, densovirus signal may derive from ingested arthropods. RNA-library detection can reflect viral transcripts or protected viral nucleic acid that remained after DNase treatment.

## Assembly and candidate selection

Non-rRNA paired reads were assembled independently for each sample with MEGAHIT. Assembly summaries were generated with QUAST, and assembled contigs were classified with MetaBuli.

Parvoviridae-classified contigs were compared with ten densovirus genomes from the previous study. Candidate selection considered nucleotide identity, query coverage, reference coverage, length, and read support.

Seven sequences were retained as near-complete current-study densovirus candidates:

- four from ABH samples;
- one from an H sample;
- two from AAH samples.

All seven candidates had 100% read breadth in their source libraries. Several contigs extend beyond the homologous region of the closest previous-study reference, so reference coverage and query coverage are reported separately.

## Read support and consensus review

Source reads were mapped back to every retained candidate. The workflow records:

- mapped read segments;
- breadth of coverage;
- mean depth;
- base and mapping quality;
- depth percentiles and terminal-to-internal depth ratios.

Potential ORF junctions were reviewed using local read depth, mismatches, insertions, and deletions. A one-nucleotide deletion candidate in `Current_Nr11` was evaluated separately. Four of five high-quality reads supported the corrected allele at the inspected junction, and the corrected 4,904-nt consensus was retained for downstream final comparisons. Because local depth is low, this correction and the affected ORF remain explicitly documented for manual review.

## Orientation and reference comparison

All study sequences were standardized to the orientation of RefSeq genome NC_031450.1, Parus major densovirus isolate PmDNV-JL. Orientation was determined by BLASTn strand support and validated after reorientation.

NC_031450.1 is used as an annotated reference and orientation anchor, not as proof that the study sequences share the same host range or species assignment.

## Protein-marker phylogeny

The phylogenetic panel contains:

- 48 external reference genomes;
- 10 genomes from Popov et al. (2025);
- 7 current-study genomes.

Reference coding sequences were extracted from GenBank annotations. Proteinortho and explicit GenBank feature information were used to curate NSP1 and SP1 marker candidates. Current- and previous-study marker proteins were reconstructed by reference-guided protein-to-genome searches.

NSP1 and SP1 were aligned separately with MAFFT. Columns with less than 70% occupancy were removed with trimAl. The two marker alignments were concatenated while retaining separate model partitions.

Three datasets were analysed:

1. **primary** — strict marker completeness criteria;
2. **study-inclusive** — retains all study genomes with sufficient marker evidence;
3. **sensitivity** — includes additional shorter or less complete marker candidates.

IQ-TREE used ModelFinder and 1,000 ultrafast bootstrap plus 1,000 SH-aLRT replicates. The selected study-inclusive partition models were:

```text
NSP1: Q.PFAM+F+R4
SP1:  Q.PFAM+F+I+R4
```

The study-inclusive topology is used for the principal manuscript figure. Primary and sensitivity trees are retained as robustness analyses. Topology comparison showed that most apparent differences occurred at weakly supported branches; the well-supported study-sequence backbone was stable.

## Whole-genome similarity

The final 18-sequence comparison comprises:

- RefSeq NC_031450.1;
- ten previous-study genomes;
- seven current-study genomes, including the reviewed `Current_Nr11` consensus.

All-versus-all BLASTn summaries report both:

- aligned nucleotide identity;
- reciprocal coverage, defined conservatively from query and reference coverage.

These metrics must be interpreted together. High identity over a short aligned region is not equivalent to high whole-genome similarity. Heatmaps therefore display identity and reciprocal coverage separately and annotate study origin.

## Conservative feature annotation

The annotation workflow searches for the five coding features represented in NC_031450.1:

- ORF5, non-structural protein;
- ORF3, non-structural protein and NSP1 candidate;
- ORF4, non-structural protein;
- ORF1, structural protein and SP1 candidate;
- ORF2, structural protein and SP2 candidate.

Only CDS calls with coherent coordinates, translation, and protein-level validation are written to provisional annotation files. A feature is omitted when the available sequence or read evidence does not support a defensible intact CDS.

For the seven current-study genomes:

- 31 CDS features passed protein validation;
- three genomes contained all five expected validated ORFs;
- four genomes retained one unresolved ORF each and were annotated conservatively.

The same conservative procedure was applied to the ten previous-study genomes. Provisional annotations are suitable for manuscript figures and pre-submission review, but they are not treated as final GenBank submissions.

## Genome-organization figure

The genome map displays validated CDS arrows, genome position, study origin, and unresolved ORFs. Text labels are omitted from the RefSeq arrows because the colour legend already identifies each ORF. The figure is descriptive and does not imply experimentally validated protein function.

## GenBank submission status

Submission has intentionally been deferred until after completion of the thesis manuscript. Before deposition, the following remain necessary:

1. final manual review of every unresolved ORF and the corrected `Current_Nr11` region;
2. confirmation of isolate names, dates, host, location, and isolation source;
3. validation of feature tables with NCBI submission tools;
4. assignment of RNA-sequencing accessions and linkage to the BioProject;
5. replacement of provisional figure labels with accession numbers after release.

The repository preserves the code needed to regenerate the evidence tables and provisional annotation set before this final review.
