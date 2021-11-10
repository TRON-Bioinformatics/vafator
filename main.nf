#!/usr/bin/env nextflow


nextflow.enable.dsl = 2

include { VAFATOR; MULTIALLELIC_FILTER } from './nf_modules/vafator'


params.help= false
params.input_files = false
params.output = "output"
params.mapping_quality = false
params.base_call_quality = false
params.skip_multiallelic_filter = false

if (params.help) {
    log.info params.help_message
    exit 0
}

// checks required inputs
if (params.input_files) {
  Channel
    .fromPath(params.input_files)
    .splitCsv(header: ['name', 'vcf', 'tumor_bams', 'normal_bams'], sep: "\t")
    .map{ row-> tuple(row.name, file(row.vcf), row.tumor_bams, row.normal_bams) }
    .set { input_files }
} else {
  exit 1, "Input file not specified!"
}

workflow {
    if (params.skip_multiallelic_filter) {
        VAFATOR(
            input_files)
    }
    else {
        MULTIALLELIC_FILTER(
            VAFATOR(
                input_files))
    }
}
