#! /usr/bin/env nextflow

/* 
 * NF RNAseq pipeline:
 * https://github.com/seqeralabs/nextflow-tutorial
 * modified by Pavel Jedlicka [jedlicka@ibp.cz] for reads from multiple samples mapping using STAR aligner and FQ quality check
 * organism: Nannochloropsis gaditana
 */

/* 
 * pipeline input parameters are called from
 * `nextflow.config` file
 */
 
log.info """\
         R N A S E Q - N F   P I P E L I N E    
         ===================================
         genome       : ${params.genome}
         annotation   : ${params.gtf}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

 
/* 
 * define the `star_index` process that create a binary `Genome` index 
 * for given genome fasta
 */

process start_index {
    cpus 2

    input:
    path genome from params.genome
    path gtf from params.gtf

    output:
    path 'star_index' into star_index_ch

    script:       
    """
    ~/Documents/Soft/STAR/source/STAR --runThreadN ${task.cpus} \
        --runMode genomeGenerate \
        --genomeDir star_index \
        --genomeSAindexNbases 11 \
        --limitGenomeGenerateRAM 124544990592 \
        --genomeFastaFiles ${genome} \
        --sjdbGTFfile ${gtf} \
        --sjdbOverhang 99    
    """
}

Channel 
    .fromFilePairs( params.reads, checkIfExists:true )
    .into { read_pairs_ch; read_pairs2_ch } 

/*
 * Run STAR to perform the quantification of expression using
 * the index and the matched read files
 */

process quantification {
    publishDir params.outdir, mode:'copy'
    cpus 2

    input:
    path index from star_index_ch
    tuple val(pair_id), path(reads) from read_pairs_ch
 
    output:
    path("${pair_id}_ReadsPerGene.out.tab") into quant_ch
 
    script:
    """
    /home/pavel/Documents/Soft/STAR/source/STAR \
        --genomeDir ${index} \
        --runThreadN ${task.cpus} \
        --limitBAMsortRAM 1239463334 \
        --readFilesCommand zcat \
        --quantMode GeneCounts \
        --readFilesIn ${reads[0]},${reads[1]} \
        --outFileNamePrefix ${pair_id}_ \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard    
    """
}

/*
 * Run fastQC to check quality of reads files
 */


process fastqc {
    tag "FASTQC on $sample_id"
    // read_pairs_ch.view()
    input:
    tuple val(sample_id), path(reads) from read_pairs2_ch

    output:
    path("fastqc_${sample_id}_logs") into fastqc_ch

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """  
}

/*
 * Create a report using multiQC for the quantification
 * and fastqc processes
 */

process multiqc {
    publishDir params.outdir, mode:'copy'

    input:
    path('*') from quant_ch.mix(fastqc_ch).collect()

    output:
    path('multiqc_report.html')

    script:
    """
    multiqc .
    """
}

workflow.onComplete { 
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
