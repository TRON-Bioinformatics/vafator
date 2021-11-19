import pysam
from cyvcf2 import VCF, Writer, Variant
import os
import vafator
import datetime
import json
from vafator.pileups import get_variant_pileup, build_variant, get_metrics


class Annotator(object):

    vafator_header = {
        "name": "vafator",
        "version": vafator.VERSION,
        "date": datetime.datetime.now().ctime(),
        "timestamp": datetime.datetime.now().timestamp(),
    }

    def __init__(self, input_vcf: str, output_vcf: str, input_bams: dict, mapping_qual_thr=0, base_call_qual_thr=29):

        self.mapping_quality_threshold = mapping_qual_thr
        self.base_call_quality_threshold = base_call_qual_thr
        self.vcf = VCF(input_vcf)
        # sets a line in the header with the command used to annotate the file
        self.vafator_header["input_vcf"] = input_vcf
        self.vafator_header["output_vcf"] = output_vcf
        self.vafator_header["bams"] = ";".join(["{}:{}".format(s, ",".join(b)) for s, b in input_bams.items()])
        self.vafator_header["mapping_quality_threshold"] = mapping_qual_thr
        self.vafator_header["base_call_quality_threshold"] = base_call_qual_thr
        self.vcf.add_to_header("##vafator_command_line={}".format(json.dumps(self.vafator_header)))
        # adds to the header all the names of the annotations
        for a in Annotator._get_headers(input_bams):
            self.vcf.add_info_to_header(a)
        self.vcf_writer = Writer(output_vcf, self.vcf)
        self.bam_readers = {s : [pysam.AlignmentFile(b, "rb") for b in bams] for s, bams in input_bams.items()}

    @staticmethod
    def _get_headers(input_bams: dict):
        headers = []

        for s, bams in input_bams.items():
            headers.append({
                'ID': "{}_af".format(s),
                'Description': "Allele frequency for the alternate alleles in the {} sample/s".format(s),
                'Type': 'Float',
                'Number': 'A'
            })
            headers.append({
                'ID': "{}_dp".format(s),
                'Description': "Total depth of coverage in the {} sample/s (independent of alleles)".format(s),
                'Type': 'Float',
                'Number': '1'
            })
            headers.append({
                'ID': "{}_ac".format(s),
                'Description': "Allele count for the alternate alleles in the {} sample/s".format(s),
                'Type': 'Integer',
                'Number': 'A'
            })

            if len(bams) > 1:
                for i, bam in enumerate(bams, start=1):
                    n = os.path.basename(bam).split(".")[0]
                    headers = headers + [
                        {'ID': "{}_af_{}".format(s, i),
                         'Description': "Allele frequency for the alternate alleles in the {} sample {}".format(s, n),
                         'Type': 'Float', 'Number': 'A'},
                        {'ID': "{}_dp_{}".format(s, i),
                         'Description': "Depth of coverage in the {} sample {} (independent of alleles)".format(s, n),
                         'Type': 'Float', 'Number': '1'},
                        {'ID': "{}_ac_{}".format(s, i),
                         'Description': "Allele count for the alternate alleles in the {} sample {}".format(s, n),
                         'Type': 'Integer', 'Number': 'A'}
                    ]
        return headers

    def _write_batch(self, batch):
        for v in batch:
            self.vcf_writer.write_record(v)

    def _add_stats(self, variant: Variant):

        global_dp = 0
        global_ac = {}
        vafator_variant = build_variant(variant)
        for sample, bams in self.bam_readers.items():
            for i, bam in enumerate(bams):
                pileups = get_variant_pileup(
                    variant=vafator_variant, bam=bam,
                    min_base_quality=self.base_call_quality_threshold,
                    min_mapping_quality=self.mapping_quality_threshold)
                coverage_metrics = get_metrics(variant=vafator_variant, pileups=pileups)
                if coverage_metrics is not None:
                    if len(bams) > 1:
                        variant.INFO["{}_af_{}".format(sample, i + 1)] = ",".join(
                            [str(self._calculate_af(coverage_metrics.ac[alt], coverage_metrics.dp)) for alt in variant.ALT])
                        variant.INFO["{}_ac_{}".format(sample, i + 1)] = ",".join([str(coverage_metrics.ac[alt]) for alt in variant.ALT])
                        variant.INFO["{}_dp_{}".format(sample, i + 1)] = coverage_metrics.dp
                    for alt in variant.ALT:
                        global_ac[alt] = global_ac.get(alt, 0) + coverage_metrics.ac[alt]
                    global_dp += coverage_metrics.dp

            variant.INFO["{}_af".format(sample)] = ",".join([str(self._calculate_af(global_ac.get(alt, 0), global_dp)) for alt in variant.ALT])
            variant.INFO["{}_ac".format(sample)] = ",".join([str(global_ac.get(alt, 0)) for alt in variant.ALT])
            variant.INFO["{}_dp".format(sample)] = global_dp

    def _calculate_af(self, ac, dp):
        return float(ac) / dp if dp > 0 else 0.0

    def run(self):
        batch = []
        variant: Variant
        for variant in self.vcf:
            # gets the counts of all bases across all BAMs
            self._add_stats(variant)

            batch.append(variant)
            if len(batch) >= 1000:
                self._write_batch(batch)
                batch = []
        if len(batch) > 0:
            self._write_batch(batch)

        self.vcf_writer.close()
        self.vcf.close()
        for _, bams in self.bam_readers.items():
            for bam in bams:
                bam.close()
