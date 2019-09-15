#!/usr/bin/env nextflow

version='0.0.1'
timestamp='20190904'


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
	System.out.prinln("	---centrifuge_db [folder]	Folder containing centrifuge index for mapping")
    exit 1
}


//Set default parameters
params.fastq = "$baseDir/*.f*q*"
params.human = "$baseDir/data/hg38.fa"
params.outdir = "results"
params.mode = "complete"
params.centrifuge_db = null

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
 genome:	${params.human}
 fastq:		${params.fastq}
 outdir:	${params.outdir}
 mode:		${params.mode}
 centrifuge_db: ${params.centrifuge_db}

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


human = file(params.human)

//Collect files
Channel
	.fromPath( params.fastq, type: 'file', checkIfExists: true)
	.into { fastqs_set; fastqs_set2; count_fastq; list_fastq}

count_fastq
	.count()
	.subscribe { println "$it file(s) detected:" }

list_fastq
	.toSortedList()
	.subscribe { println it }

fastqs_set
	.merge(fastqs_set2)
	.set {fastqs}

//validate fastq

process validate {
conda 'bioconda::fqtools=2.0'
cpus 1
memory '4 GB'

input:
set file(input),val(filename) from fastqs

output:
set file(input), val(filename), stdout into validation_results

"""
fqtools validate $input
"""

}


passed_validation = Channel.create()
failed_validation = Channel.create()


validation_results.choice(passed_validation, failed_validation) {
    validation_results -> validation_results[2] == "OK\n" ? 0 : 1
}

passed_validation
	.into { passed_fastqs; passed_validation2; passed_fastqs2 }

failed_validation
	.subscribe { failed_validation -> println "The following sample failed fastq validation: " + failed_validation[1] }

passed_validation2	
        .ifEmpty { error "No samples passed validation. Exiting now." }	
	.subscribe { passed_validation2 -> println "The following sample passed fastq validation: " + passed_validation2[1] }
	
        


//fastqc


process fastqc {
conda 'bioconda::fastqc=0.11.8'
//container 'biocontainers/fastqc:v0.11.5_cv4'
cpus 1
memory '2 GB'

input:
set file(fastq),val(path),val(passed) from passed_fastqs

output:
file '*_fastqc.zip' into fastqc

"""
fastqc $fastq
"""

}



//multiqc pre preprocess

process multiqc_fastqc {
conda 'bioconda::multiqc=1.7'
//container 'ewels/multiqc'
publishDir params.outdir

input:
file '*_fastqc.zip' from fastqc.collect()
output:
file '*.html'

"""
multiqc *_fastqc.zip
"""

}

// Trim nanopore adapters. Porechop seems to be deprecated...still need to trim adapters without barcode? qcat currently doesn't support this.
process demultiplex {
	conda 'bioconda::qcat'
	cpus 8
	memory '4 GB'
	echo true
	
	input:
	set file(input),val(path),val(passed) from passed_fastqs2	

	output:
	set val(input.baseName), file('output/*.fastq') into chopped mode flatten
	
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
	set val(filename), file(input) from chopped

	output:
	set val(filename), val(input.baseName), stdout into high_complexity

	"""
	prinseq-lite.pl -min_qual_mean 7 -lc_method dust -lc_threshold 7 -min_len 250 -fastq $input -out_good stdout -out_bad null
	
	"""


}


//Build centrifuge index

process build_centrifuge_index {
	conda 'bioconda::centrifuge bioconda::blast'
	memory '300 GB'
	cpus 32
	publishDir './Centrifuge_index'
	echo true 

	output:
	file('*.cf') into centrifuge_index
        
	when:
	params.centrifuge_db == null

        
	"""
	centrifuge-download -P 32 -o taxonomy taxonomy
	centrifuge-download -P 32 -o library -m -a "Any" -d "archaea,viral,fungi,protozoa" refseq > seqid2taxid.map
	centrifuge-download -P 32 -o library -m -d "bacteria" refseq >> seqid2taxid.map
	centrifuge-download -P 32 -o library -m -d "vertebrate_mammalian" -a "Chromosome" -t 9606 -c 'reference genome' refseq >> seqid2taxid.map
	centrifuge-download -P 32 -o library -m contaminants >> seqid2taxid.map
	cat library/*/*.fna > input-sequences.fna
	centrifuge-build -p 32 --conversion-table seqid2taxid.map --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp input-sequences.fna refseq_microbial
        """

}



/*
//Centrifuge!
//The awk filters by read coverage. $6 is hit length and $7 is query length
process centrifuge {
	conda 'bioconda::centrifuge'
	memory '100 GB'
	cpus 32

	input:
	set val(filename), val(adapter), stdin from high_complexity
	
	if params.centrifuge_db == null {
	file('*.cf') from centrifuge_index
	}
	else {
	file(params.centrifuge_index)
	}


	output:
	set val(filename), val(adapter), file('centrifuge.class') into centrifuge_results		

	"""
	centrifuge -p 32 -U - -k 1 --min-hitlen 16 --host-taxids 9606 -x $reference | awk -F "\t" 'BEGIN {OFS = FS} NR==1{print} NR>1{if(($6/$7>0.9)&&($6>0)) print}' > centrifuge.class
	"""


}
*/

/*
//Recentrifuge for visualization, subtracting controls
process recentrifuge_getdb {
	conda 'bioconda::recentrifuge'
	cpus 1

	output:
	file('*') into recentrifuge_db


	
	"""
	retaxdump
	"""
}

//I suppose could also use Krona but can filter with recentrifuge
process recentrifuge {
	conda 'bioconda::recentrifuge'
	cpus 8

	input:
	set val(filename), val(adapter), file(input) from centrifuge_results.collect()
	file('*') from recentrifuge_db
	output:

	"""
	rcf -f $input -x 9606 -e TSV -s NORMA
	"""

}

//Extract reads of species detected for coverage calculation

//Assemble species with high coverage


//AMR

//Species specific pipelines
//TB resistance
//Toxin detection
//MLST

process mlst {
	conda '/conda_configs/mlst.yaml'
	cpus 1
	



*/
//Plot and report



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
