#!/usr/bin/env nextflow

version='0.0.1'
timestamp='20190904'
params.version = null
params.help = null


//Gets starting time		
sysdate = new java.util.Date()

/**
	Prints version when asked for
*/
if (params.version) {
	System.out.println("")
	System.out.println("BUGSEQ PIPELINE - Version: $version ($timestamp)")
	exit 1
}

/**
	Prints help when asked for
*/

if (params.help) {
	System.out.println("")
	System.out.println("BUGSEQ PIPELINE - Version: $version ($timestamp)")
	System.out.println("Usage: ")
	System.out.println("   nextflow run bugseq.nf --fastqs [Fastq files] --outdir path  ")
	System.out.println("")
	System.out.println("Mandatory arguments:")
	System.out.println("    --fastq   *.f*q      File path(s) of fastq file(s). Can be uncompressed or gzipped.")
	System.out.println("Optional arguments:")

    exit 1
}


//Set default parameters
params.fastq = "$baseDir/*.f*q*"
params.human = "$baseDir/data/hg38.fa"
params.outdir = "results"
params.mode = "complete"
params.centrifuge_db = null
params.control_fastq = null




//Checking user-defined parameters	
if (params.mode != "complete" && params.mode != "filter" && params.mode != "process") {
	exit 1, "Mode not available. Choose any of <filter, process, complete>"
}	


//Collect system information
osname = System.getProperty("os.name") //Operating system name
osarch = System.getProperty("os.arch") //Operating system architecture
osversion = System.getProperty("os.version") //Operating system version
osuser = System.getProperty("user.name") //User's account name

javaversion = System.getProperty("java.version") //Java Runtime Environment version
javaVMname = System.getProperty("java.vm.name") //Java Virtual Machine implementation name
javaVMVersion = System.getProperty("java.vm.version") //Java Virtual Machine implementation version


log.info """\
 B U G S E Q - N F   P I P E L I N E
 ===================================
 fastqs:	${params.fastq}
 outdir:	${params.outdir}
 mode:		${params.mode}
 centrifuge_db: ${params.centrifuge_db}
 control_fastq:	${params.control_fastq}
 
++++++++++++++++++++++++++++++++++++

 Start time: ${sysdate}
 Operating System:
	name:         $osname
	architecture: $osarch
	version:      $osversion
 Java
	version: $javaversion
	Java Virtual Machine: $javaVMname ; version: $javaVMVersion
 Nextflow:
	version:   $nextflow.version	
	build:     $nextflow.build	
	timestamp: $nextflow.timestamp
 ===================================
 """
.stripIndent()

// GET READY TO RUMBLE (PREP REFS)


//This is to calculate metagenome size
process get_genome_sizes {
	cpus 1
	memory '500 MB'
	
	output:
	file('species_genome_size.txt') into ncbi_genome_size

	script:

	"""
	wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/species_genome_size.txt.gz -O - | gunzip > species_genome_size.txt
	"""
}


//This is for the metagenomic assembly taxonomic classification
Channel
	.fromPath(params.centrifuge_db + '.*.cf', type: 'file')
	.ifEmpty { error "Please specify the centrifuge index. Basename up to .1.cf should be included." }
	.collect()
	.set { centrifuge_index }
	



//Get recentrifuge taxonomy
process recentrifuge_getdb {
	conda '/home/schorlton/.conda/envs/recentrifuge/' //Would need to upload this package
	cpus 1

	output:
	file('taxdump/') into recentrifuge_db


	
	"""
	retaxdump
	"""
}



// OK Let's actually start processing some files now


//Collect files
Channel
	.fromPath( params.fastq, type: 'file', checkIfExists: true)
	.tap { count_fastq; list_fastq}
	.map { tuple(it,it,"case") }
	.set { case_fastqs }

count_fastq
	.count()
	.subscribe { println "$it case file(s) detected:" }

list_fastq
	.toSortedList()
	.subscribe { println it }


if ( params.control_fastq ) {
Channel
	.fromPath( params.control_fastq, type: 'file', checkIfExists: true)
	.tap { count_fastq_control; list_fastq_control}
	.map { tuple(it,it,"control") }
	.set { control_fastqs }

control_count=count_fastq_control
	.count()
	.subscribe { println "$it control file(s) detected:" }

list_fastq_control
	.toSortedList()
	.subscribe { println it }

case_fastqs
	.mix(control_fastqs)
	.set { fastqs }


}
else {
fastqs = case_fastqs
control_count=0
}




//validate fastq

process validate {
conda 'bioconda::fqtools=2.0'
cpus 1
memory '4 GB'

input:
set file(input),val(filename),val(case_control) from fastqs

output:
set file(input), val(filename), val(case_control), stdout into validation_results

"""
fqtools validate $input
"""

}


passed_validation = Channel.create()
failed_validation = Channel.create()


validation_results.choice(passed_validation, failed_validation) {
    validation_results -> validation_results[3] == "OK\n" ? 0 : 1
}

passed_validation
	.tap { passed_fastqs; passed_fastqs2 }
        .ifEmpty { exit 1,  "No samples passed validation. Exiting now." }	
	.subscribe { passed_validation -> println "The following samples passed fastq validation: " + passed_validation[1] }

failed_validation
	.subscribe { failed_validation -> exit 1, println "The following sample failed fastq validation: " + failed_validation[1] }	
	


//fastqc


process fastqc {
conda 'bioconda::fastqc=0.11.8'
//container 'biocontainers/fastqc:v0.11.5_cv4'
cpus 1
memory '2 GB'

input:
set file(fastq),val(path),val(case_control), val(passed) from passed_fastqs

output:
file '*_fastqc.zip' into fastqc

"""
fastqc $fastq
"""

}




// Trim nanopore adapters. Porechop seems to be deprecated...still need to trim adapters without barcode? qcat currently doesn't support this.
process demultiplex {
	conda 'bioconda::qcat future'
	cpus 8
	memory '4 GB'
	echo true
	
	input:
	set file(input),val(path),val(case_control),val(passed) from passed_fastqs2	

	output:
	set val(path), file('output/*.fastq'), val(case_control) into chopped mode flatten
	
	"""
	qcat -f $input --trim -t 8 -b output
	"""

}



//Based on CrumpIT and EPI2ME cutoffs.
process filter_low_complexity {
	conda 'bioconda::prinseq'
	cpus 2
	memory '2 GB'
	//publishDir 'results/filtered/'
	
	input:
	set val(filename), file(input), val(case_control) from chopped

	output:
	set val(filename), val(input.baseName), val(case_control), file('out.fastq') into high_complexity
	set val(filename), file('out.fastq') into high_complexity2

	"""
	prinseq-lite.pl -min_qual_mean 7 -lc_method dust -lc_threshold 7 -min_len 250 -fastq $input -out_good out -out_bad null
	
	"""


}

//Detect experiment type here. Is it amplicon, isolate or metagenomic?




//If metagenomic:
//Taxonomic read level classification!

process read_classify {
	conda 'bioconda::minimap2'
	memory '125 GB'
	cpus 32
	publishDir params.outdir

	input:
	set val(filename), val(adapter), val(case_control), file(input) from high_complexity
	file minimap2_index
	
	output:
	set val(filename), val(adapter), val(case_control), stdout into centrifuge_results

	script:
	index_base = centrifuge_index[0].toString() - ~/.\d.cfl?/
	
	"""
	minimap2 -t 32 -x $index_base $input

	"""


}

// Refine alignments

process refine_alignments {
	conda 'bioconda::pathoscope'


}


//I suppose could also use Krona but can filter contaminants with recentrifuge


cases_recentrifuge = Channel.create()
controls_recentrifuge = Channel.create()

centrifuge_results_withCov
	.choice(cases_centrifuge,controls_centrifuge) {
	centrifuge_results_withCov -> centrifuge_results_withCov[2] == "case" ? 0 : 1
	}

cases_recentrifuge	
	.into { cases_recentrifuge_files ; cases_recentrifuge_command ; cases_recentrifuge_tofilter }

controls_recentrifuge
	.into { controls_recentrifuge_files; controls_recentrifuge_command }



//Need to decide whether to change -m (mintaxa) which is the reads before collapsing.
process recentrifuge {
	conda '/home/schorlton/.conda/envs/recentrifuge/'
	cpus 8
	publishDir params.outdir	

	input:
	
	file cases from cases_recentrifuge_files.collect { it[3] }
	val case_cmd from cases_recentrifuge_command.collect{[' -g ', it[3].baseName+'.class'].flatten()}
	
	file controls from controls_recentrifuge_files.collect { it[3] }
	val control_cmd from controls_recentrifuge_command.collect{[' -g ', it[3].baseName+'.class'].flatten()}
	
	val control_count	

	file('*') from recentrifuge_db
	

	output:
	
	file('*.class.rcf.data.tsv') into recentrifuge_out
	
	script:	
	
	if (control_count>0){	
	"""	
	rcf --format "TYP:TSV, TID:3, LEN:7, SCO:9, UNC:0" -x 9606 -s GENERIC -e TSV -y 30 -c $control_count${control_cmd.join()}${case_cmd.join()}
	"""
	}
	else {
	"""
	rcf --format "TYP:TSV, TID:3, LEN:7, SCO:9, UNC:0" -x 9606 -s GENERIC -e TSV -y 30${case_cmd.join()}
	"""
}
}


/*


//Filter human reads and contamination from the recentrifuge results, calculate assembly size, assemble, then run centrifuge on assemblies and match up assemblies with centrifuge results. Assemblies will be needed for later steps.

process metagenomic_assembly {
	conda 'bioconda::flye'
	cpus 32
	memory { fastq.size() < 1.GB ? 100.GB : 350.GB }
	
	publishDir params.outdir
	//echo true
	errorStrategy { task.attempt < 4 ? 'retry' : 'ignore' }
	maxRetries 3
	input:
	set val(case_filename), val(case_adapter), file(fastq), val(size) from to_assemble

	output:
	set val(case_filename), val(case_adapter), file('flye_out/assembly.fasta') into metagenomes_quast, metagenomes2, metagenomes_amr, metagenomes_toxins
	
	script:
	mysize=size.trim().toInteger()/task.attempt
	"""
	flye --meta --plasmid --nano-raw $fastq --genome-size ${mysize.round(0)}B -o flye_out --threads 32
	"""

}


//Assess assembly size with quast/metaquast
//eventually I suppose should use reference genomes extracted from centrifuge taxids. Currently metaquast uses 16s alignment to figure out reference genomes...

process metaquast {

	conda 'bioconda::quast'
	cpus 8
	memory '8 GB'

	publishDir params.outdir	
	input:
	set val(filename), val(adapter), file("${filename.baseName}_${adapter}.fasta") from metagenomes_quast
	output:
	file('quast_results/latest/report.tsv') into metagenomes_quast_results
	script:
	"""
	metaquast --max-ref-number 0 --threads 8 ${filename.baseName}_${adapter}.fasta
	"""
	

}


//Centrifuge the metagenomic assembly
process metagenome_centrifuge {
	conda 'bioconda::centrifuge'
	memory '200 GB'
	cpus 32
	

	input:
	set val(filename), val(adapter), file(input) from metagenomes2
	file centrifuge_index
	
	output:
	set val(filename), val(adapter), stdout into centrifuge_metagenome_results

	script:
	index_base = centrifuge_index[0].toString() - ~/.\d.cfl?/
	
	"""
	centrifuge -p 32 -f -U $input -k 1 -x $index_base

	"""


}

//Calculate contig coverage
//The awk filters by read coverage. $6 is hit length and $7 is query length
//Can't get centrifuge to throw error when piping to awk if the reference is bad
//32630 is contaminant taxid in centrifuge as per https://github.com/DaehwanKimLab/centrifuge/blob/6cc874e890249f6d8b479bd41b41e131f12c6796/centrifuge-download#L259

//at some point it would also be nice to remove contamination pre-assembly based on recentrifuge results

process calculate_metagenomeCov {
	cpus 1
	memory '1 GB'
	publishDir params.outdir
	input:
	set val(filename), val(adapter), stdin from centrifuge_metagenome_results

	output:
	set val(filename), val(adapter), file('*.class') into centrifuge_metagenome_withCov
		
	
	"""
 	awk -F"\t" '{(NR > 1) ? \$(NF+1)=\$6/\$7*100 : \$(NF+1)="readCov"}1' OFS="\t" > ${filename.baseName}_${adapter}_meta.class

	"""

}

//alternatives include abricate and rgi. Phenotype info in RGI not great. next step: combine species info for SNP detection
process search_assembly_amr {
	conda 'bioconda::staramr'
	cpus 32
	publishDir params.outdir
	input:
	set val(case_filename), val(case_adapter), file(assembly) from metagenomes_amr

	output:
	set val(case_filename), val(case_adapter), file('resfinder.txt') into staramar
	script:

	"""
	staramr search $assembly --pid-threshold 90 --percent-length-overlap-resfinder 90 --nprocs 32 --output-summary /dev/null --output-resfinder resfinder.txt
	"""	

}

//abricate for toxins

toxin_fasta = Channel.fromPath(workflow.projectDir + '/refs/toxins.fasta')

process search_toxins {
	conda '/conda_configs/abricate.yaml'
	cpus 8
	memory '4 GB'
	publishDir params.outdir

	input:
	set val(case_filename), val(case_adapter), file(assembly) from metagenomes_toxins 
	output:
	set val(case_filename), val(case_adapter), file('out.tab') into abricate_toxins

	script:
	"""
	abricate --threads 8 --mincov 90 --db toxins --datadir refs/ $assembly > out.tab
	"""


}
*/

/*


//Species specific pipelines

process mlst {
	conda '/conda_configs/mlst.yaml'
	cpus 1


}
	
//TB resistance

//MLST and cg/wgMLST if possible
//B cereus vs anthracis
//Salmonella serotyping: SISTR
//H flu serotyping: hicap
//E coli serotyping
//Strep pneumo serotyping
//resistance for strep pneumo andn gonorrheae: https://www.biorxiv.org/content/10.1101/403204v2.full


*/
//Plot and report

/*
//multiqc pre preprocess

process multiqc_fastqc {
conda 'bioconda::multiqc=1.7'
//container 'ewels/multiqc'
publishDir params.outdir + '/QC/'

input:
file '*_fastqc.zip' from fastqc.collect()
file '*' from metagenomes_quast_results.collect()
output:
file '*.html'

"""
multiqc -n "BUGSEQ_QC" *
"""

}


process report {

input:
set amr_genepresence

output:
file('report.html')
file('report.pdf')

script:




}
*/
workflow.onComplete {
    println """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}
