#!/usr/bin/env nextflow

params.help= false
params.input_files = false
params.output = false
params.mapping_quality = false
params.base_call_quality = false

def helpMessage() {
    log.info"""
Usage:
    nextflow run main.nf --input_files input_files

This workflow implements the annotation of a tumor/normal pair VCF with the
variant allele frequencies, counts and depth of coverage extracted from the
corresponding BAM files.

Input:
    * input_files: the path to a tab-separated values file containing in each
    row the sample name, path to the VCF file, a comma-separated list of normal
    BAMs and a comma-separated list of normal BAMs
    The input file does not have header!
    Example input file:
    identifier	/path/to/your/file.vcf  /path/to/your/normal_replica1.bam,/path/to/your/normal_replica2.bam /path/to/your/tumor_replica1.bam,/path/to/your/tumor_replica2.bam

Optional input:
    * output: the folder where to publish output

Output:
    * Annotated VCF file
    """
}

if (params.help) {
    helpMessage()
    exit 0
}

publish_dir = "output"
if (params.output) {
  publish_dir = params.output
}

// checks required inputs
if (params.input_files) {
  Channel
    .fromPath(params.input_files)
    .splitCsv(header: ['name', 'vcf', 'normal_bams', 'tumor_bams'], sep: "\t")
    .map{ row-> tuple(row.name, file(row.vcf), row.normal_bams, row.tumor_bams) }
    .set { input_files }
} else {
  exit 1, "Input file not specified!"
}

process vafator {
    cpus 1
    memory '4g'
    //container 'biocontainers/bcftools'
    module 'anaconda/3/2019'
    tag "${name}"
    publishDir "${publish_dir}/${name}", mode: "copy"

    input:
    	set name, file(vcf), normal_bams, tumor_bams from input_files

    output:
      set val("${name}"), val("${publish_dir}/${name}/${vcf.baseName}.vaf.vcf") into output_files
      file("${vcf.baseName}.vaf.vcf") into annotated_vcf

    script:
    normal_bams_param = normal_bams?.trim() ? "--normal-bams " + normal_bams.split(",").join(" ") : ""
    tumor_bams_param = tumor_bams?.trim() ? "--tumor-bams " + tumor_bams.split(",").join(" ") : ""
    mapping_quality_param = params.mapping_quality ? "--mapping-quality " + params.mapping_quality : ""
    base_call_quality_param = params.base_call_quality ? "--base_call-quality " + params.base_call_quality : ""
    """
    #source /home/priesgof/src/vafator/venv/bin/activate
    /home/priesgof/src/vafator/venv/bin/vafator --input-vcf ${vcf} --output-vcf ${vcf.baseName}.vaf.vcf ${normal_bams_param} ${tumor_bams_param} ${mapping_quality_param} ${base_call_quality_param}
    """
}

output_files
	.map {it.join("\t")}
	.collectFile(name: "${publish_dir}/annotated_vcfs.txt", newLine: true)
