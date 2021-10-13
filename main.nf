#!/usr/bin/env nextflow


nextflow.enable.dsl = 2

include { VAFATOR; MULTIALLELIC_FILTER } from './nf_modules/vafator'


params.help= false
params.input_files = false
params.output = "output"
params.mapping_quality = false
params.base_call_quality = false
params.skip_multiallelic_filter = false
params.prefix = false

if (params.help) {
    log.info params.help_message
    exit 0
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
