# Vafator

Annotate a VCF file with allele frequency, allele counts and depth of coverage for alternate alleles.

Outputs a VCF with the following annotations for tumor and normal:
```
chr1    12345       .       A       G       .       PASS  tumor_af=0.0;tumor_ac=0;tumor_dp=89;normal_af=0.0196078431372549;normal_ac=1;normal_dp=51
```

Both tumor and normal BAMs are optional, it can only annotate with the tumor BAM or viceversa.
If more than one BAM is provided for either the tumor or the normal the frequencies and counts are added up.