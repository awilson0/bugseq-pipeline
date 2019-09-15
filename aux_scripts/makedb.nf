#!/usr/bin/env nextflow
params.centrifuge_db = null
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

