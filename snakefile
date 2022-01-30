##########################################################################
# This is the Snakefile of a workflow to annotate plant genomes using    #
# InterProScan and BLAST (functional annotation).                        #
#                                                                        #
#    1) Extract proteins                                                 #
#    2) Remove Stop-Codons                                               #
#    3) Identify Genes with internal Stop                                #
#    4) InterProScan                                                     #
#    5) BLAST                                                            #
#                                                                        #
##########################################################################

##########################################################################
########################## SNAKEMAKE PARAMETERS ##########################
##########################################################################

configfile: "config.yaml"

import os

REF_PATH=os.path.dirname(config["REF"]) + "/"
REF_GFF="/augustus.hints.gff3"

OUT_DIR=config["ANN_DIR"] + "/" + config["REF_NAME"] + "/"

PROT_DIR=os.path.dirname(config["PROT"]) + "/"
PROT_FASTA=os.path.basename(config["PROT"])
PROT_NAME, PROT_EXT=os.path.splitext(PROT_FASTA) #EXAMPLE: protein.faa --> PROT_NAME = protein, PROT_EXT = faa

##########################################################################
############################ SNAKEMAKE RULES #############################
##########################################################################

rule all:
    input:
        OUT_DIR + "augustus.noInternalStop.gff.stats",
        OUT_DIR + "augustus.interpro.blast.gff.stats",
        OUT_DIR + "eggnog/augustus.interpro.blast.eggnog.emapper.gff"

### Functional annotation using InterproScan and Blast ###

rule extract_proteins:
    """Extract protein sequences from reference assembly based on the augustus gtf file, with the transcript ID as fasta header"""
    input:
        ref = config["REF"],
        gff = config["BRA_DIR"] + REF_GFF
    output:
        fasta = OUT_DIR + "augustus.protein.fasta"
    params:
        dir = OUT_DIR
    shell:
        """
        cd {params.dir}
        module load cufflinks
        gffread -y {output.fasta} -O -E -L -F -g {input.ref} {input.gff}
        """


rule remove_stop_codons_end:
    """Remove '.' from the end of lines that are used as stop codons in gffread, which causes errors when running blast and interproscan to fail"""
    input:
        rules.extract_proteins.output.fasta
    output:
        fasta = OUT_DIR + "augustus.protein.noStopEnd.fasta"
    params:
        dir = OUT_DIR
    shell:
        """
        cd {params.dir}
        sed 's/\.*$//g' {input} > {output.fasta}
        """

rule identify_genes_with_internal_stop:
    """Make a list of genes that contain internal stop codons"""
    input:
        rules.remove_stop_codons_end.output.fasta
    output:
        OUT_DIR + "augustus.protein.genesInternalStop.list"
    params:
        dir = OUT_DIR,
        code_dir = "scripts"
    shell:
        """
        module unload python3
        module load biopython
        python3 {params.code_dir}/identify_genes_internal_stop_codons.py {input} {output}
        """

rule remove_internal_stop_augustus_gff:
    """Remove genes with internal stop codons (coded as '.') from the augustus gff file"""
    input:
        genelist = rules.identify_genes_with_internal_stop.output,
        gff = config["BRA_DIR"] + REF_GFF
    output:
        fixed_gff = OUT_DIR + "augustus.noInternalStop.gff"
    params:
        dir = OUT_DIR
    shell:
        """
        awk '$0="ID="$0' {input.genelist} | grep -v -f - {input.gff} > intermediate.gff
        echo  "##gff-version 3" > {output.fixed_gff} && cat intermediate.gff >> {output.fixed_gff}
        rm intermediate.gff
        """

rule no_internal_stop_gff_stats:
    """Get some simple stats from the gff file"""
    input:
        rules.remove_internal_stop_augustus_gff.output.fixed_gff
    output:
        gff_stats = OUT_DIR + "augustus.noInternalStop.gff.stats"
    params:
        dir = OUT_DIR
    shell:
        """
        cd {params.dir}
        module load genometools
        grep -v '# gffread' {input} | gt stat - > {output.gff_stats}
        """

rule remove_internal_stop_protein_fasta:
    """Remove genes with internal stop codons (coded as '.') from the protein fasta file"""
    input:
        genelist = rules.identify_genes_with_internal_stop.output,
        fasta = rules.remove_stop_codons_end.output.fasta
    output:
        tmp = temp(OUT_DIR + "augustus.protein.genesInternalStop.IDs.list"),
        fasta = OUT_DIR + "augustus.protein.noInternalStop.fasta"
    params:
        dir = OUT_DIR,
        code_dir = "scripts"
    shell:
        """
        awk '$0=$0"_t"' {input.genelist} > {output.tmp}
        module unload python3
        module load biopython
        python3 {params.code_dir}/remove_fasta_entries.py {input.fasta} {output.tmp} {output.fasta}
        """

rule interproscan:
    """Search for pathways, families, domains, sites, repeats, structural domains and other sequence features in the InterproScan database"""
    input:
        rules.remove_internal_stop_protein_fasta.output.fasta 
    output:
        tsv_ipr = OUT_DIR + "augustus.protein.noInternalStop.fasta.tsv",
        gff3_ipr = OUT_DIR + "augustus.protein.noInternalStop.fasta.gff3",
        xml_ipr = OUT_DIR + "augustus.protein.noInternalStop.fasta.xml"
    params:
        dir = OUT_DIR
    log: 
    	"interproscan.log"
    shell:
        """
        cd {params.dir}
        /localmirror/monthly/interpro/interproscan-5.52-86.0/interproscan.sh -i {input} -t p -dp -pa --goterms --iprlookup --pathways --cpu 16 &> {log}
        """

rule blast_db:
    """Generate BLAST database"""
    input: config["PROT"]
    output:
        phr = config["PROT"] + ".phr",
        pin = config["PROT"] + ".pin",
        psq = config["PROT"] + ".psq"
    params:
        dir = PROT_DIR
    shell:
        """
        cd {params.dir}
        module load ncbiblastplus
        makeblastdb -in {input} -dbtype prot
        """

rule blast:
    """BLASTp of proteins from augustus against uniprot database"""
    input:
        fasta = rules.remove_internal_stop_protein_fasta.output.fasta,
        blast_db_idx = rules.blast_db.output,
        blast_db = config["PROT"] 
    output:
        blast = OUT_DIR + "augustus.noInternalStop.blastp.out"
    params:
        dir = OUT_DIR
    threads: 8
    shell:
        """
        cd {params.dir}
        module load ncbiblastplus
        blastp -num_threads {threads} -db {input.blast_db} -query {input.fasta} -evalue 10e-5 -outfmt 6 -out {output.blast}
        """

rule update_gff_blast:
    """Add BLAST results to gff file"""
    input:
        gff = config["BRA_DIR"] + REF_GFF,
        ipr = rules.interproscan.output.tsv_ipr,
        blast = rules.blast.output.blast,
        blast_db = config["PROT"]
    output:
        gff = OUT_DIR + "augustus.interpro.blast.gff"
    params:
        dir = OUT_DIR,
        agatdir = OUT_DIR + "augustus_interpro_blast"
    conda:
    	"envs/agat.yaml"
    shell:
        """
	module unload perl
        cd {params.dir}
        agat_sp_manage_functional_annotation.pl -f {input.gff} -b {input.blast} --db {input.blast_db} -i {input.ipr} --output {params.agatdir} &&
        cp {params.agatdir}/augustus.hints.gff3 {output.gff}
        """

rule ipr_blast_gff_stats:
    """Get some simple stats from the gff file incl. InterProScan and/or BLAST results"""
    input:
        rules.update_gff_blast.output.gff
    output:
        blast_stats = OUT_DIR + "augustus.interpro.blast.gff.stats"
    params:
        dir = OUT_DIR
    shell:
        """
        cd {params.dir}
        module load genometools
        grep -v '# gffread' {input} | gt stat - > {output.blast_stats}
        """
        
rule emapper:
    """Add KEGG Annotation results for Pathway analysis using eggnogmapper"""
    input:
        fa = rules.remove_internal_stop_protein_fasta.output.fasta,
        gff = rules.update_gff_blast.output.gff
    output:
        enogg = OUT_DIR + "eggnog/augustus.interpro.blast.eggnog.emapper.gff"
    params:
        out = "augustus.interpro.blast.eggnog",
        dir = OUT_DIR + "eggnog",
	name = config["REF_NAME"]
    shell:
        """
        cd {params.dir}
        module load eggnogmapper
        emapper.py -m diamond -i {input.fa} --tax_scope 33090 --target_taxa 33090 --report_orthologs --decorate_gff {input.gff} --output {params.out} --override --cpu 0 &&
        cd ..
        cp eggnog/augustus.interpro.blast.eggnog.emapper.gff .
	mv augustus.interpro.blast.eggnog.emapper.gff {params.name}_annotated.gff
        """
        


