# Annotation batch size for VCF writing
BATCH_SIZE = 10000

# IUPAC ambiguity codes treated as ambiguous bases in pileup metrics
AMBIGUOUS_BASES = ['N', 'M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B']

# Header templates: (suffix, description, type, number)
# {sample} is substituted at generation time
_HEADER_TEMPLATES = [
    ("af",       "Allele frequency for the alternate alleles in the {sample} sample/s",                                                                              "Float",   "A"),
    ("dp",       "Total depth of coverage in the {sample} sample/s (independent of alleles)",                                                                        "Float",   "1"),
    ("ac",       "Allele count for the alternate alleles in the {sample} sample/s",                                                                                  "Integer", "A"),
    ("n",        "Allele count for ambiguous bases (any IUPAC ambiguity code is counted) in the {sample} sample/s",                                                  "Integer", "1"),
    ("pu",       "Probability of an undetected mutation given the observed supporting reads (AC), the observed total coverage (DP) and the provided tumor purity in the {sample} sample/s", "Float", "A"),
    ("pw",       "Power to detect a somatic mutation as described in Absolute given the observed total coverage (DP) and the provided tumor purity and ploidies in the {sample} sample/s",  "Float", "1"),
    ("k",        "Minimum number of supporting reads, k, such that the probability of observing k or more non-reference reads due to sequencing error is less than the defined FPR in the {sample} sample/s", "Float", "1"),
    ("eaf",      "Expected VAF considering the purity and ploidy/copy number in the {sample} sample/s",                                                              "Float",   "1"),
    ("bq",       "Median base call quality of the reads supporting each allele in the {sample} sample/s",                                                            "Float",   "R"),
    ("mq",       "Median mapping quality of the reads supporting each allele in the {sample} sample/s",                                                              "Float",   "R"),
    ("pos",      "Median position within the read of the reads supporting each allele in the {sample} sample/s",                                                     "Float",   "R"),
    ("rsmq",     "Rank sum test comparing the MQ distributions supporting the reference and the alternate in the {sample} sample/s",                                 "Float",   "A"),
    ("rsmq_pv",  "Rank sum test p-value for MQ distributions in the {sample} sample/s. The null hypothesis is that there is no difference between the distributions", "Float",  "A"),
    ("rsbq",     "Rank sum test comparing the BQ distributions supporting the reference and the alternate in the {sample} sample/s",                                 "Float",   "A"),
    ("rsbq_pv",  "Rank sum test p-value for BQ distributions in the {sample} sample/s. The null hypothesis is that there is no difference between the distributions", "Float",  "A"),
    ("rspos",    "Rank sum test comparing the position distributions supporting the reference and the alternate in the {sample} sample/s",                            "Float",   "A"),
    ("rspos_pv", "Rank sum test p-value for position distributions in the {sample} sample/s. The null hypothesis is that there is no difference between the distributions", "Float", "A"),
]

# eaf is not produced per-replicate
_REPLICATE_HEADER_TEMPLATES = [t for t in _HEADER_TEMPLATES if t[0] != "eaf"]