from cyvcf2 import Variant
from pysam.libcalignmentfile import IteratorColumnRegion


def get_insertion_metrics(variant: Variant, pileups: IteratorColumnRegion):
    ac = 0
    dp = 0
    position = variant.POS
    insertion_length = len(variant.ALT[0]) - len(variant.REF)
    try:
        pileup = next(pileups)
        for r in pileup.pileups:
            dp += 1
            if r.indel > 0:
                # read with an insertion
                start = r.alignment.reference_start
                for cigar_type, cigar_length in r.alignment.cigartuples:
                    if cigar_type in [0, 2, 3, 7, 8]:  # consumes reference M, D, N, =, X
                        start += cigar_length
                        if start > position:
                            break
                    elif cigar_type == 1:  # does not count I
                        if start == position and cigar_length == insertion_length:
                            ac += 1
    except StopIteration:
        # no reads
        pass
    return ac, dp


def get_deletion_metrics(variant: Variant, pileups: IteratorColumnRegion):
    ac = 0
    dp = 0
    position = variant.POS
    deletion_length = len(variant.REF) - len(variant.ALT[0])
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
                            ac += 1
                        else:
                            start += cigar_length
                    if start > position:
                        break
    except StopIteration:
        # no reads
        pass
    return ac, dp


def _initialize_empty_count(bases_counts, base):
    if bases_counts.get(base) is None:
        bases_counts[base] = 0
    return bases_counts


def get_snv_metrics(variant: Variant, pileups: IteratorColumnRegion):
    dp = 0
    ac = {alt: 0 for alt in variant.ALT}
    try:
        pileup = next(pileups)
        for s in pileup.get_query_sequences():
            dp += 1
            if s in ac:
                ac[s] = ac[s] + 1
    except StopIteration:
        # no reads
        pass
    return ac, dp
