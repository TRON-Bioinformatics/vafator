params.cpus = 1
params.memory = "4g"
params.output = ""
params.mapping_quality = false
params.base_call_quality = false
params.skip_multiallelic_filter = false
params.prefix = false


process VAFATOR {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy"

    input:
    tuple val(name), file(vcf), val(normal_bams), val(tumor_bams)

    output:
    tuple val("${name}"), file("${vcf.baseName}.vaf.vcf"), emit: annotated_vcf

    script:
    normal_bams_param = normal_bams?.trim() ? "--normal-bams " + normal_bams.split(",").join(" ") : ""
    tumor_bams_param = tumor_bams?.trim() ? "--tumor-bams " + tumor_bams.split(",").join(" ") : ""
    mq_param = params.mapping_quality ? "--mapping-quality " + params.mapping_quality : ""
    bq_param = params.base_call_quality ? "--base-call-quality " + params.base_call_quality : ""
    prefix_param = params.prefix ? "--prefix " + params.prefix : ""
    """
    vafator \
    --input-vcf ${vcf} \
    --output-vcf ${vcf.baseName}.vaf.vcf ${normal_bams_param} ${tumor_bams_param} ${mq_param} ${bq_param} ${prefix_param}
    """
}


process MULTIALLELIC_FILTER {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy"

    input:
    tuple val(name), file(vcf)

    output:
    file("${vcf.baseName}.filtered_multiallelics.vcf") //, emit: filtered_vcf

    script:
    """
    multiallelics-filter --input-vcf ${vcf} --output-vcf ${vcf.baseName}.filtered_multiallelics.vcf
    """
}
