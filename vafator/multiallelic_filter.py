from cyvcf2 import VCF, Writer
import vafator
import datetime
import json
import random


class MultiallelicFilter(object):

    vafator_header = {
        "name": "multiallelic-filter",
        "version": vafator.VERSION,
        "date": datetime.datetime.now().ctime(),
        "timestamp": datetime.datetime.now().timestamp(),
    }

    multiallelic_annotation_name = "multiallelic"

    def __init__(self, input_vcf, output_vcf, tumor_sample_name='tumor'):
        """
        :param input_vcf: the input VCF file
        :param output_vcf: the file path of the output VCF
        """
        self.tumor_sample_name = tumor_sample_name
        self.vcf = VCF(input_vcf)
        # sets a line in the header with the command used to annotate the file
        self.vafator_header["input_vcf"] = input_vcf
        self.vafator_header["output_vcf"] = output_vcf
        self.vcf.add_to_header("##vafator_command_line={}".format(json.dumps(self.vafator_header)))
        # adds to the header all the names of the annotations
        self.vcf.add_info_to_header({'ID': self.multiallelic_annotation_name, 'Type': 'String', 'Number': '.',
                 'Description': "Indicates multiallelic variants filtered and their frequencies if any (e.g.: T,0.12)"})
        self.vcf_writer = Writer(output_vcf, self.vcf)

    def run(self):
        batch = []
        prev_variant = None
        for variant in self.vcf:
            if prev_variant is None:
                # first time stores the variant and moves on
                prev_variant = variant
                continue
            # considers only SNVs with same chromosome, position and reference
            if variant.is_snp and variant.CHROM == prev_variant.CHROM and variant.POS == prev_variant.POS and \
                    variant.REF == prev_variant.REF:
                af1 = self.get_tumor_af(prev_variant)
                af2 = self.get_tumor_af(variant)
                # keeps the variant with the highest AF
                if af1 < af2:
                    self.set_multiallelic_annotation(variant, prev_variant.ALT[0], af1)
                    prev_variant = variant
                elif af1 > af2:
                    self.set_multiallelic_annotation(prev_variant, variant.ALT[0], af2)
                elif af1 == af2:
                    # chooses a variant at random
                    alt1 = variant.ALT[0]
                    alt2 = prev_variant.ALT[0]
                    prev_variant = random.sample([variant, prev_variant], k=1)[0]
                    self.set_multiallelic_annotation(prev_variant, alt1 if alt1 != prev_variant.ALT[0] else alt2, af1)
                continue

            # write previous variant and stores current for next iteration
            batch.append(prev_variant)
            prev_variant = variant

            if len(batch) >= 1000:
                self._write_batch(batch)
                batch = []
        # do not forget to add the last variant!
        batch.append(prev_variant)
        self._write_batch(batch)

        self.vcf_writer.close()
        self.vcf.close()

    def set_multiallelic_annotation(self, variant, alt, af):
        variant.INFO[self.multiallelic_annotation_name] = \
            ",".join(variant.INFO.get(self.multiallelic_annotation_name, "").split(",") + [alt, str(af)])

    def get_tumor_af(self, prev_variant):
        return prev_variant.INFO.get("{}_af".format(self.tumor_sample_name), 0.0)

    def _write_batch(self, batch):
        for v in batch:
            self.vcf_writer.write_record(v)
