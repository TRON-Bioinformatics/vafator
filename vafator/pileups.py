from dataclasses import dataclass
from typing import List
from cyvcf2 import Variant
from pysam.libcalignmentfile import IteratorColumnRegion, AlignmentFile


@dataclass
class VafatorVariant:
    chromosome: str
    position: int
    reference: str
    alternative: List[str]

    def is_snp(self):
        return len(self.reference) == 1 and len(self.alternative[0]) == 1

    def is_insertion(self):
        return len(self.reference) == 1 and len(self.alternative[0]) > 1

    def is_deletion(self):
        return len(self.alternative[0]) == 1 and len(self.reference) > 1


def build_variant(variant: Variant) -> VafatorVariant:
    return VafatorVariant(
        chromosome=variant.CHROM,
        position=variant.POS,
        reference=variant.REF,
        alternative=variant.ALT
    )


def get_variant_pileup(
        variant: VafatorVariant, bam: AlignmentFile, min_base_quality, min_mapping_quality) -> IteratorColumnRegion:
    chromosome = variant.chromosome
    position = variant.position
    # this function returns the pileups at all positions covered by reads covered the queried position
    # approximately +- read size bp
    return bam.pileup(contig=chromosome, start=position - 1, stop=position, truncate=True,
                      min_base_quality=min_base_quality,
                      min_mapping_quality=min_mapping_quality)


@dataclass
class CoverageMetrics:
    ac: dict
    dp: int


def get_metrics(variant: VafatorVariant, pileups: IteratorColumnRegion) -> CoverageMetrics:
    if variant.is_snp():
        return get_snv_metrics(variant, pileups)
    elif variant.is_insertion():
        return get_insertion_metrics(variant, pileups)
    elif variant.is_deletion():
        return get_deletion_metrics(variant, pileups)
    return None


def get_insertion_metrics(variant: VafatorVariant, pileups: IteratorColumnRegion) -> CoverageMetrics:
    ac = {alt.upper(): 0 for alt in variant.alternative}
    dp = 0
    position = variant.position
    insertion_length = len(variant.alternative[0]) - len(variant.reference)
    insertion = variant.alternative[0][1:]
    try:
        pileup = next(pileups)
        for r in pileup.pileups:
            dp += 1
            if r.indel > 0:
                # read with an insertion
                start = r.alignment.reference_start
                relative_position = 0
                for cigar_type, cigar_length in r.alignment.cigartuples:
                    if cigar_type in [0, 2, 3, 7, 8]:  # consumes reference M, D, N, =, X
                        start += cigar_length
                        if start > position:
                            break
                    if cigar_type in [0, 1, 4, 7, 8]:  # consumes query M, I, S, =, X
                        relative_position += cigar_length
                    if cigar_type == 1:  # does not count I
                        insertion_in_query = r.alignment.query[relative_position:relative_position + insertion_length]
                        if start == position and cigar_length == insertion_length and \
                                insertion == insertion_in_query:
                            ac[variant.alternative[0].upper()] = ac[variant.alternative[0].upper()] + 1
    except StopIteration:
        # no reads
        pass
    return CoverageMetrics(ac=ac, dp=dp)


def get_deletion_metrics(variant: VafatorVariant, pileups: IteratorColumnRegion) -> CoverageMetrics:
    ac = {alt.upper(): 0 for alt in variant.alternative}
    dp = 0
    position = variant.position
    deletion_length = len(variant.reference) - len(variant.alternative[0])
    try:
        pileup = next(pileups)
        for r in pileup.pileups:
            dp += 1
            if r.indel < 0:
                # read with an deletion
                start = r.alignment.reference_start
                for cigar_type, cigar_length in r.alignment.cigartuples:
                    if cigar_type in [0, 3, 7, 8]:  # consumes reference M, N, =, X
                        start += cigar_length
                    elif cigar_type == 2:  # D
                        if start == position and cigar_length == deletion_length:
                            ac[variant.alternative[0].upper()] = ac[variant.alternative[0].upper()] + 1
                        else:
                            start += cigar_length
                    if start > position:
                        break
    except StopIteration:
        # no reads
        pass
    return CoverageMetrics(ac=ac, dp=dp)


def _initialize_empty_count(bases_counts, base):
    if bases_counts.get(base) is None:
        bases_counts[base] = 0
    return bases_counts


def get_snv_metrics(variant: VafatorVariant, pileups: IteratorColumnRegion) -> CoverageMetrics:
    dp = 0
    ac = {alt.upper(): 0 for alt in variant.alternative}
    ac[variant.reference.upper()] = 0
    try:
        pileup = next(pileups)
        for s in pileup.get_query_sequences():
            dp += 1
            base_upper_case = s.upper()
            if base_upper_case in ac:
                ac[base_upper_case] = ac[base_upper_case] + 1
    except StopIteration:
        # no reads
        pass
    return CoverageMetrics(ac=ac, dp=dp)
