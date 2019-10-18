# bugseq-pipeline
A pipeline to analyze clinical microbiology nanopore sequencing data

## Quick start
```
git clone https://github.com/schorlton/bugseq-pipeline.git
cd bugseq-pipeline
nextflow main.nf --fastq in.fq --outdir output_dir
```
## Outline
* [Requirements](#requirements)
* [Input](#input)
* [Output](#output)

## Requirements
* [Nextflow](https://www.nextflow.io)
* [Conda](https://docs.conda.io/en/latest/miniconda.html)
* A powerful [computer](https://en.wikipedia.org/wiki/Computer) (eg. 128GB RAM or more)

## Input
A nanopore basecalled fastq. Can be any library version (R7-10). Can be barcoded or not. Can be amplicon data (eg. 16S/ITS), isolate data (eg. a colony of *Staphylococcus aureus*), or clinical metagenomic data (eg. the sputum of a patient).

## Output
An interactive html file and a static pdf summary file.
Example [here](https://en.wikipedia.org/wiki/Special:Random)

## All options
```
usage: nextflow main.nf --fastq 'in.fq' --outdir 'out_dir' [options...]
options:
  # Input options
  --fastq
  --control
  
  # Output options
  --outdir
  
  # Pipeline options
  --isolate
  --metagenome
  --16S
  --ITS
  --skipQC
  --skipTyping
  --skipAMR
```

## Pipeline overview
1. User inputs basecalled nanopore fastq reads
2. BugSeq validates the fastq file (fqtools) and determines if it's truly nanopore data
3. Next, the fastq undergoes quality assessment with FastQC and results combined with multiqc
4. Reads are adapter trimmed and demultiplexed (poretools)
5. Trimmed and demultiplexed reads undergo experiment type detection to determine if this is amplicon data (eg. 16S/ITS), cultured isolate data or metagenomic data.

### Isolate data
1. Genome assembly (flye)
2. Taxonomic classification of assembly

### Metagenomic data
1. Read-level taxonomic classification to species level
2. Correction for control samples to identify significant pathogens in the cases only
3. Metagenome assembly
4. Taxonomic binning of species within metagenome

### Pathogen specific analyses
  * Phenotypic antimicrobial resistance prediction
  * MLST (when public scheme avaiable)
  * cgMLST/wgMLST (when public scheme available)
  * Serotyping (when applicable)
  * Other old-school typing (eg. Spoligotyping for *M. tuberculosis*)
  * Toxin detection (when clinically relevant)
  * Phylogenetic tree building (when multiple isolates inputted)
  
### Changelog
0.0.1
