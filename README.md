# bugseq-pipeline
BugSeq automatically analyzes clinical microbiology nanopore sequencing data from start to finish. This includes taxonomic classification of reads, antimicrobial resistance prediction and detailed subtyping for public health purposes. It was created during [hackseq 2019](https://www.hackseq.com/hackseq19#descriptions)

## Rationale
Modern clinical microbiology techniques take a day to grow and identify an organism, and another day to determine antimicrobial susceptibilities. Yet, clinical trials show that patients with septic shock have a ~7% increase in mortality for every hour delay in appropriate antimicrobial therapy. Furthermore, patients with rare or novel infections may never have the etiology of their illness diagnosed as traditional techniques can not pick up their identify their infection. Metagenomic nanopore sequencing has the potential to drastically speed up the diagnosis and characterization of infections, potentially including novel pathogens, enabling better patient outcomes. Recovering pathogen genomes provides a vast amount of clinically useful information, such as whether a patient's *E. coli* is susceptible to ceftriaxone, whether a patient's *V. cholerae* is toxigenic, or whether the *M. tuberculosis* between two patients are likely to be epidemiologically linked.


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
An interactive html file and a static pdf summary file. These files will show a taxonomic classification of the percentage makeup of organisms in each patient sample. This could include viral, bacterial, fungal and protozoal organisms. For each of these organisms, the presence of antimicrobial resistance genes will be visualized, along with the predicted phenotype for the antimicrobial drug associated with these genes. 

Example [here](https://en.wikipedia.org/wiki/Special:Random)

## All options
```
usage: nextflow main.nf --fastq 'in.fq' --outdir 'out_dir' [options...]
options:
  # Input options
  --fastq
  --control_fastq            Control sample fastqs for complex subtraction from cases. Can use regex patterns or specify multiple file separated by commas.
  
  # Output options
  --outdir
  
  # Pipeline options
  --isolate                 Input file(s) are from isolate sequencing
  --metagenome              Input file(s) are metagenomic samples
  --16S                     Input file(s) are 16S amplicon sequencing data
  --ITS                     Input file(s) are ITS amplicon sequencing data
  --skipQC                  Skip quality assessment and read trimming
  --skipTyping              Skip public health typing analysis
  --skipAMR                 Skip AMR prediction step
  --meanQ [7]               Reads with mean quality below this value will be filtered from analysis
  --minLength [250]         Reads with length below this threshold will be filtered from analysis
```

## Pipeline overview
1. User inputs basecalled nanopore fastq reads
2. BugSeq validates the fastq file (fqtools) and determines if it's truly nanopore data
3. Next, the fastq undergoes quality assessment with FastQC and results combined with multiqc
4. Reads are adapter trimmed and demultiplexed (qcat)
5. Reads are quality and length filtered
6. Trimmed and demultiplexed reads undergo experiment type detection to determine if this is amplicon data (eg. 16S/ITS), cultured isolate data or metagenomic data (magic..., including sourmash)

### Isolate data
1. Genome assembly (Flye)
2. Taxonomic classification of assembly (minimap2 + Pathoscope ID)

### Metagenomic data
1. Read-level taxonomic classification to species level
2. Correction for control samples to identify significant pathogens in the cases only
3. Metagenome assembly (metaFlye)
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
