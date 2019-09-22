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
	System.out.println("	---centrifuge_db [folder]	Folder containing centrifuge index for mapping")
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


//Get resistance database. CARD database could be an issue with licensing but let's run with it for now.

process get_card {

	cpus 1
	memory '1 GB'

	output:
	file('nucleotide_fasta_protein_homolog_model.fasta') into card_fasta

	"""
	wget https://card.mcmaster.ca/latest/data -O - | tar -xj
	"""
}


toxin_fasta = Channel.fromPath(workflow.projectDir + '/refs/toxins.fasta')

// Build minimap2 resistance and toxin database
process build_resistance_db {
	conda 'bioconda::minimap2'
	cpus 8
	memory '16 GB'
	
	input:
	file card_fasta
	file toxin_fasta
	output:
	file('nuc_prot_homolog.mmi') into resistance_db
	file('toxins.mmi') into toxin_db
	script:
	"""
	minimap2 -t 8 -x map-ont -d nuc_prot_homolog.mmi $card_fasta
	minimap2 -t 8 -x map-ont -d toxins.mmi $toxin_fasta
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
	set val(filename), val(input.baseName), val(case_control), file('out.fastq') into high_complexity, high_complexity_AST
	set val(filename), file('out.fastq') into high_complexity2

	"""
	prinseq-lite.pl -min_qual_mean 7 -lc_method dust -lc_threshold 7 -min_len 250 -fastq $input -out_good out -out_bad null
	
	"""


}




//Centrifuge!
//Still need to incorporate nextflow generated db if none specified
process centrifuge {
	conda 'bioconda::centrifuge'
	memory '200 GB'
	cpus 32
	publishDir params.outdir

	input:
	set val(filename), val(adapter), val(case_control), file(input) from high_complexity
	file centrifuge_index
	
	output:
	set val(filename), val(adapter), val(case_control), stdout into centrifuge_results

	script:
	index_base = centrifuge_index[0].toString() - ~/.\d.cfl?/
	
	"""
	centrifuge -p 32 -U $input -k 1 --min-hitlen 16 --host-taxids 9606 -x $index_base

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
		
	
	"""
 	awk -F"\t" '{(NR > 1) ? \$(NF+1)=\$6/\$7*100 : \$(NF+1)="readCov"}1' OFS="\t" > ${filename.baseName}_${adapter}.class

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
	.into { cases_centrifuge_files ; cases_centrifuge_command ; cases_centrifuge_tofilter }

controls_centrifuge
	.into { controls_centrifuge_files; controls_centrifuge_command }



//Need to decide whether to change -m (mintaxa) which is the reads before collapsing.
process recentrifuge {
	conda '/home/schorlton/.conda/envs/recentrifuge/'
	cpus 8
	publishDir params.outdir	

	input:
	
	file cases from cases_centrifuge_files.collect { it[3] }
	val case_cmd from cases_centrifuge_command.collect{[' -g ', it[3].baseName+'.class'].flatten()}
	
	file controls from controls_centrifuge_files.collect { it[3] }
	val control_cmd from controls_centrifuge_command.collect{[' -g ', it[3].baseName+'.class'].flatten()}
	
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




toextract = cases_centrifuge_tofilter.join(high_complexity2)


//For now, only remove human and contaminant reads before assembly. Maybe at a later date, also remove contamination reads pre-assembly. Will be ?harder to calculate assembly size though.
//Only require 10% read coverage for human and contaminant reads
//Also estimate metagenome size here
process extract_foreground_reads {

	conda 'bioconda::seqtk'
	cpus 1
	memory '1 GB'
	publishDir params.outdir

	input:
	file ncbi_genome_size
	set val(case_filename), val(case_adapter), val(case_control_case), file(centrifuge_out), file(fastq_in) from toextract


	output:
	set val(case_filename), val(case_adapter), file('keepers.fastq'), stdout into to_assemble
	file('*')
	script:
	
	"""	
	awk -F"\t" '(NR>3 && \$9>40 && \$3!=9606 && \$3!=32630){print \$3}' $centrifuge_out | sort -u > mytaxids.txt

	awk -F"\t" '(NR>3 && !(\$9>10 && ((\$3=9606) || (\$3=32630)))){print \$1}' $centrifuge_out > toextract_reads.txt
	seqtk subseq $fastq_in toextract_reads.txt > keepers.fastq

	awk 'NR == FNR{a[\$1] = \$4;next}; {sum+=a[\$1]} END {print sum}' $ncbi_genome_size mytaxids.txt


	"""

}


//Filter human reads and contamination from the recentrifuge results, calculate assembly size, assemble, then run centrifuge on assemblies and match up assemblies with centrifuge results. Assemblies will be needed for later steps.

process metagenomic_assembly {
	conda 'bioconda::flye'
	cpus 32
	memory { fastq.size() < 1.GB ? 100.GB : 300.GB }
	
	//publishDir params.outdir
	//echo true
	errorStrategy { task.attempt < 4 ? 'retry' : 'ignore' }
	maxRetries 3
	input:
	set val(case_filename), val(case_adapter), file(fastq), val(size) from to_assemble

	output:
	set val(case_filename), val(case_adapter), file('flye_out/assembly.fasta') into metagenomes
	
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

	input:
	set val(filename), val(adapter), file("${filename.baseName}_${adapter}.fasta") from metagenomes
	output:
	file('quast_results/latest/report.tsv') into metagenomes_quast
	script:
	"""
	metaquast --max-ref-number 0 --threads 8 ${filename.baseName}_${adapter}.fasta
	"""
	

}



//Map for AMR and toxins
// Filter anbiotic resistance gene calls by depth and breadth of coverage
process map_ast {
	conda 'bioconda::minimap2 bioconda::samtools'
	cpus 16
	memory '8 GB'
	publishDir params.outdir

	input:
	 set val(filename), val(input.baseName), val(case_control), file(reads) from high_complexity_AST
	file resistance_db
	file toxin_db
	output:
	file('*')
	
	script:

	"""
	minimap2 --sam-hit-only --secondary=no -t 14 -a $resistance_db $reads | samtools sort --threads 2 > amr.bam
	samtools depth -a amr.bam | awk -F"\t" '{seen[\$1]+=\$3; count[\$1]++; if(\$3>0)breadth[\$1]++} END{for (x in seen) if(((seen[x]/count[x])>0.8) && ((breadth[x]/count[x])>0.9))print x, seen[x]/count[x], breadth[x]/count[x]}' > found_amr.txt

	minimap2 --sam-hit-only --secondary=no -t 14 -a $toxin_db $reads | samtools sort --threads 2 > toxin.bam
	samtools depth -a toxin.bam | awk -F"\t" '{seen[\$1]+=\$3; count[\$1]++; if(\$3>0)breadth[\$1]++} END{for (x in seen) if(((seen[x]/count[x])>0.8) && ((breadth[x]/count[x])>0.9))print x, seen[x]/count[x], breadth[x]/count[x]}' > found_toxins.txt


	"""


}







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



*/
//Plot and report


//multiqc pre preprocess

process multiqc_fastqc {
conda 'bioconda::multiqc=1.7'
//container 'ewels/multiqc'
publishDir params.outdir + '/QC/'

input:
file '*_fastqc.zip' from fastqc.collect()
file '*' from metagenomes_quast.collect()
output:
file '*.html'

"""
multiqc -n "BUGSEQ_QC" *
"""

}

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
