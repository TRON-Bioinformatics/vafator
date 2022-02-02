from collections import Counter
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
    return bam.pileup(contig=chromosome, start=position - 1, stop=position,
                      truncate=True,                    # returns only this column
                      max_depth=0,                      # disables maximum depth
                      min_base_quality=min_base_quality,
                      min_mapping_quality=min_mapping_quality)


@dataclass
class CoverageMetrics:
    ac: dict
    dp: int


def get_metrics(variant: VafatorVariant, pileups: IteratorColumnRegion) -> CoverageMetrics:
    if variant.is_snp():
        return get_snv_metrics(pileups)
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
    alt_upper = variant.alternative[0].upper()
    try:
        pileups = [p for p in next(pileups).pileups if p.indel > 0]
        dp += len(pileups)
        for p in pileups:
            # read with an insertion
            start = p.alignment.reference_start
            relative_position = 0
            for cigar_type, cigar_length in p.alignment.cigartuples:
                if cigar_type in [0, 2, 3, 7, 8]:  # consumes reference M, D, N, =, X
                    start += cigar_length
                    if start > position:
                        break
                if cigar_type in [0, 1, 4, 7, 8]:  # consumes query M, I, S, =, X
                    relative_position += cigar_length
                if cigar_type == 1:  # does not count I
                    insertion_in_query = p.alignment.query[relative_position:relative_position + insertion_length]
                    if start == position and cigar_length == insertion_length and \
                            insertion == insertion_in_query:
                        ac[alt_upper] = ac[alt_upper] + 1
    except StopIteration:
        # no reads
        pass
    return CoverageMetrics(ac=ac, dp=dp)


def get_deletion_metrics(variant: VafatorVariant, pileups: IteratorColumnRegion) -> CoverageMetrics:
    ac = {alt.upper(): 0 for alt in variant.alternative}
    dp = 0
    position = variant.position
    deletion_length = len(variant.reference) - len(variant.alternative[0])
    alt_upper = variant.alternative[0].upper()
    try:
        pileups = [p for p in next(pileups).pileups if p.indel < 0]
        dp += len(pileups)
        for r in pileups:
            # read with a deletion
            start = r.alignment.reference_start
            for cigar_type, cigar_length in r.alignment.cigartuples:
                if cigar_type in [0, 3, 7, 8]:  # consumes reference M, N, =, X
                    start += cigar_length
                elif cigar_type == 2:  # D
                    if start == position and cigar_length == deletion_length:
                        ac[alt_upper] = ac[alt_upper] + 1
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


def get_snv_metrics(pileups: IteratorColumnRegion) -> CoverageMetrics:
    dp = 0
    ac = Counter()
    try:
        pileup = next(pileups)
        bases = [s.upper() for s in pileup.get_query_sequences()]
        dp += len(bases)
        ac = ac + Counter(bases)
    except StopIteration:
        # no reads
        pass
    return CoverageMetrics(ac=ac, dp=dp)
