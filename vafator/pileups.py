from collections import Counter
from dataclasses import dataclass
from typing import Union
from cyvcf2 import Variant
from pysam.libcalignmentfile import IteratorColumnRegion, AlignmentFile
from vafator.tests.utils import VafatorVariant


def is_snp(variant: Variant):
    return len(variant.REF) == 1 and len(variant.ALT[0]) == 1


def is_insertion(variant: Variant):
    return len(variant.REF) == 1 and len(variant.ALT[0]) > 1


def is_deletion(variant: Variant):
    return len(variant.ALT[0]) == 1 and len(variant.REF) > 1


def get_variant_pileup(
        variant: Union[Variant, VafatorVariant], bam: AlignmentFile, min_base_quality, min_mapping_quality) -> IteratorColumnRegion:
    position = variant.POS
    # this function returns the pileups at all positions covered by reads covered the queried position
    # approximately +- read size bp
    return bam.pileup(contig=variant.CHROM, start=position - 1, stop=position,
                      truncate=True,                    # returns only this column
                      max_depth=0,                      # disables maximum depth
                      min_base_quality=min_base_quality,
                      min_mapping_quality=min_mapping_quality)


@dataclass
class CoverageMetrics:
    ac: dict
    dp: int


def get_metrics(variant: Variant, pileups: IteratorColumnRegion) -> CoverageMetrics:
    if is_snp(variant):
        return get_snv_metrics(pileups)
    elif is_insertion(variant):
        return get_insertion_metrics(variant, pileups)
    elif is_deletion(variant):
        return get_deletion_metrics(variant, pileups)
    return None


def get_insertion_metrics(variant: Variant, pileups: IteratorColumnRegion) -> CoverageMetrics:
    ac = {alt.upper(): 0 for alt in variant.ALT}
    dp = 0
    position = variant.POS
    insertion_length = len(variant.ALT[0]) - len(variant.REF)
    insertion = variant.ALT[0][1:]
    alt_upper = variant.ALT[0].upper()
    try:
        pileups = next(pileups).pileups
        dp += len(pileups)
        for p in pileups:
            if p.indel > 0:
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


def get_deletion_metrics(variant: Variant, pileups: IteratorColumnRegion) -> CoverageMetrics:
    ac = {alt.upper(): 0 for alt in variant.ALT}
    dp = 0
    position = variant.POS
    deletion_length = len(variant.REF) - len(variant.ALT[0])
    alt_upper = variant.ALT[0].upper()
    try:
        pileups = next(pileups).pileups
        dp += len(pileups)
        for p in pileups:
            if p.indel < 0:
                # read with a deletion
                start = p.alignment.reference_start
                for cigar_type, cigar_length in p.alignment.cigartuples:
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


def get_snv_metrics(pileups: IteratorColumnRegion) -> CoverageMetrics:
    try:
        pileup = next(pileups)
        bases = [s.upper() for s in pileup.get_query_sequences()]
        dp = len(bases)
        ac = Counter(bases)
    except StopIteration:
        # no reads
        dp = 0
        ac = Counter()
    return CoverageMetrics(ac=ac, dp=dp)
