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
	System.out.prinln("	---centrifuge_db [folder]	Folder containing centrifuge index for mapping")
    exit 1
}


//Set default parameters
params.fastq = "$baseDir/*.f*q*"
params.human = "$baseDir/data/hg38.fa"
params.outdir = "results"
params.mode = "complete"
params.centrifuge_db = null
params.control_fastq = null
params.control_centrifuge = null



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
 control_centrifuge:	${params.control_centrifuge}
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
	set val(filename), val(input.baseName), val(case_control), stdout into high_complexity

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

if (params.centrifuge_db != null) {
	Channel
		.fromPath(params.centrifuge_db + '.*.cf', type: 'file')
		.ifEmpty { error "Please specify the centrifuge index. Basename up to .1.cf should be included." }
		.collect()
		.set { centrifuge_index }
	
}


//Centrifuge!
//Still need to incorporate nextflow generated db if none specified
process centrifuge {
	conda 'bioconda::centrifuge'
	memory '200 GB'
	cpus 32
	publishDir params.outdir

	input:
	set val(filename), val(adapter), val(case_control), stdin from high_complexity
	file centrifuge_index
	
	output:
	set val(filename), val(adapter), val(case_control), stdout into centrifuge_results

	script:
	index_base = centrifuge_index[0].toString() - ~/.\d.cfl?/
	
	"""
	centrifuge -p 32 -U - -k 1 --min-hitlen 16 --host-taxids 9606 -x $index_base

	"""


}

//Calculate read coverage and create a channel without human reads for assembly
//The awk filters by read coverage. $6 is hit length and $7 is query length
//Can't get centrifuge to throw error when piping to awk if the reference is bad
//32630 is contaminant taxid in centrifuge as per https://github.com/DaehwanKimLab/centrifuge/blob/6cc874e890249f6d8b479bd41b41e131f12c6796/centrifuge-download#L259

//at some point it would also be nice to remove contamination pre-assembly based on recentrifuge results

process calculate_readCov {
	cpus 1
	memory '1 GB'
	publishDir params.outdir
	input:
	set val(filename), val(adapter), val(case_control), stdin from centrifuge_results

	output:
	set val(filename), val(adapter), val(case_control), file('*.class') into centrifuge_results_withCov
	set val(filename), val(adapter), val(case_control), file('nohuman.txt') into no_human_class
	
	"""
 	awk -F"\t" '{(NR > 1) ? \$(NF+1)=\$6/\$7*100 : \$(NF+1)="readCov"}1' OFS="\t" > ${filename.baseName}_${adapter}.class
	awk -F"\t" '(NR==1){print}; (NR>1 && \$9>40 && \$3!=9606 && \$3!=32630){print}' ${filename.baseName}_${adapter}.class > nohuman.txt
	"""

}


//Recentrifuge for visualization, subtracting controls
process recentrifuge_getdb {
	conda '/home/schorlton/.conda/envs/recentrifuge/' //Would need to upload this package
	cpus 1

	output:
	file('taxdump/') into recentrifuge_db


	
	"""
	retaxdump
	"""
}

//I suppose could also use Krona but can filter contaminants with recentrifuge


cases_centrifuge = Channel.create()
controls_centrifuge = Channel.create()

centrifuge_results_withCov
	.choice(cases_centrifuge,controls_centrifuge) {
	centrifuge_results_withCov -> centrifuge_results_withCov[2] == "case" ? 0 : 1
	}

cases_centrifuge
	.into { cases_centrifuge_files ; cases_centrifuge_command }

controls_centrifuge
	.into { controls_centrifuge_files; controls_centrifuge_command }


process recentrifuge {
	conda '/home/schorlton/.conda/envs/recentrifuge/'
	cpus 8
	publishDir params.outdir	
	echo true
	input:
	set val(case_filename), val(case_adapter), val(case_control_case), file("${case_filename.baseName}_${case_adapter}.class") from cases_centrifuge_files.collect()
	val case_cmd from cases_centrifuge_command.collect{[' -g ', it[3].baseName+'.class'].flatten()}
	
	set val(control_filename), val(control_adapter), val(control_case_control), file("${control_filename.baseName}_${control_adapter}.class") from controls_centrifuge_files.collect()
	val control_cmd from controls_centrifuge_command.collect{[' -g ', it[3].baseName+'.class'].flatten()}
	
	val control_count	

	file('*') from recentrifuge_db
	

	output:
	file('*')
	
	script:	
	
	if (control_count>0){	
	"""	
	rcf --format "TYP:TSV, TID:3, LEN:7, SCO:9, UNC:0" -x 9606 -s GENERIC -e TSV -y 40 -c $control_count${control_cmd.join()}${case_cmd.join()}
	"""
	}
	else {
	"""
	rcf --format "TYP:TSV, TID:3, LEN:7, SCO:9, UNC:0" -x 9606 -s GENERIC -e TSV -y 40${case_cmd.join()}
	"""
}
}


process get_genomes_sizes {
	cpus 1
	memory '500 MB'
	
	output:
	file('species_genome_size.txt') into ncbi_genome_size

	script:

	"""
	wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/species_genome_size.txt.gz -O - | gunzip > species_genome_size.txt
	"""
}


process estimage_metagenome_size {
	echo true
	
	input:
	file ncbi_genome_size
	set val(filename), val(adapter), val(case_control), file(classification) from no_human_class
	output:
	
	val(filename), val(adapter), val(case_control), file(classification), val(stdout) into metagenome_size
	
	script:

	"""
	cut -f3 -d'	' $classification | tail -n +2 | sort -u > mytaxids.txt
	awk 'NR == FNR{a[\$1] = \$4;next}; {sum+=a[\$1]} END {print sum}' $ncbi_genome_size mytaxids.txt
	"""

}
/*
process extract_foreground_reads {

	conda 'bioconda::seqtk'
	cpus 1
	memory '1 GB'


	input:	
	val(filename), val(adapter), val(case_control), file(classification), val(size) from metagenome_size
	output:
	val(filename), val(adapter), val(case_control), file(classification), val(size), file('to_assemble.fastq') into to_assemble

	script:
	seq







}
*/
/*

//Filter human reads and contamination from the recentrifuge results, calculate assembly size, assemble, then run centrifuge on assemblies and match up assemblies with centrifuge results. Assemblies will be needed for later steps.

process metagenomic_assembly {
	conda 'bioconda::flye'
	cpus 32
	memory '300 GB'

	output:

	"""
	flye --meta --plasmid --nano-raw $input --genome-size $metagenome_size -o flye_out --threads 32
	"""



}


//Assess assembly size with metaquast


//Get ?NCBI database. CARD database could be an issue with licensing

//Build minimap2 reference

//Map

//Calculate coverage


//AMR

//Species specific pipelines

process mlst {
	conda '/conda_configs/mlst.yaml'
	cpus 1


}
	
//TB resistance
//Toxin detection: diphtheria, other corynes for diphtheria toxin, Clostridium botulinum for botulinum toxin, Vibrio
//MLST and cg/wgMLST if possible
//B cereus vs anthracis
//Salmonella serotyping: SISTR
//H flu serotyping: hicap
//E coli serotyping
//Strep pneumo serotyping



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
