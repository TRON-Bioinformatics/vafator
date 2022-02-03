from cyvcf2 import VCF


def _get_count_variants(input_file):
    vcf = VCF(input_file)
    n_variants = 0
    for v in vcf:
        n_variants += 1
    vcf.close()
    return n_variants


def _get_info_fields(input_file):
    vcf = VCF(input_file)
    return [h.info().get("ID") for h in vcf.header_iter() if h['HeaderType'] == 'INFO']


class VafatorVariant:

    def __init__(self, chromosome, position, reference, alternative):
        self.CHROM = chromosome
        self.POS = position
        self.REF = reference
        self.ALT = alternative
