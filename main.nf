#!/usr/bin/env nextflow


if (params.help) {
    log.info params.help_message
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
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${publish_dir}/${name}", mode: "copy"

    input:
    	set name, file(vcf), normal_bams, tumor_bams from input_files

    output:
      set val("${name}"), file("${vcf.baseName}.vaf.vcf") into annotated_vcf

    script:
    normal_bams_param = normal_bams?.trim() ? "--normal-bams " + normal_bams.split(",").join(" ") : ""
    tumor_bams_param = tumor_bams?.trim() ? "--tumor-bams " + tumor_bams.split(",").join(" ") : ""
    mapping_quality_param = params.mapping_quality ? "--mapping-quality " + params.mapping_quality : ""
    base_call_quality_param = params.base_call_quality ? "--base-call-quality " + params.base_call_quality : ""
    """
    vafator --input-vcf ${vcf} --output-vcf ${vcf.baseName}.vaf.vcf ${normal_bams_param} ${tumor_bams_param} ${mapping_quality_param} ${base_call_quality_param}
    """
}

if (!params.skip_multiallelic_filter) {
    process multiallelic_filter {
        cpus params.cpus
        memory params.memory
        tag "${name}"
        publishDir "${publish_dir}/${name}", mode: "copy"

        input:
            set name, file(vcf) from annotated_vcf

        output:
          set val("${name}"), val("${publish_dir}/${name}/${vcf.baseName}.filtered_multiallelics.vcf") into output_files
          file("${vcf.baseName}.filtered_multiallelics.vcf") into filtered_vcf

        script:
        """
        multiallelics-filter --input-vcf ${vcf} --output-vcf ${vcf.baseName}.filtered_multiallelics.vcf
        """
    }
}
