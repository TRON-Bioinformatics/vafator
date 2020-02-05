from cyvcf2 import VCF


def _get_count_variants(input_file):
    vcf = VCF(input_file)
    n_variants = 0
    for v in vcf:
        n_variants += 1
    vcf.close()
    return n_variants