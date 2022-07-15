from cyvcf2 import VCF


def _get_count_variants(input_file):
    vcf = VCF(input_file)
    n_variants = 0
    for v in vcf:
        n_variants += 1
    vcf.close()
    return n_variants


def _get_mutation_at_position(input_file, chromosome, position):
    variant = None
    vcf = VCF(input_file)
    for v in vcf(): # we cannot query by specific positions as this requires a tabix index
        if v.CHROM == chromosome and v.POS == position:
            variant = v
            break
    vcf.close()
    return variant


def _get_info_fields(input_file):
    vcf = VCF(input_file)
    return [h.info().get("ID") for h in vcf.header_iter() if h['HeaderType'] == 'INFO']


def _get_annotation_values(input_file, annotation):
    vcf = VCF(input_file)
    values = []
    for v in vcf:
        values.append(v.INFO.get(annotation))
    return values


class VafatorVariant:

    def __init__(self, chromosome, position, reference, alternative):
        self.CHROM = chromosome
        self.POS = position
        self.REF = reference
        self.ALT = alternative
