import os.path
import unittest

from unassigner.ani import RefseqAssembly

DATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")
SUMMARY_FP = os.path.join(DATA_DIR, "assembly_summary.txt")


class RefseqAssemblyTests(unittest.TestCase):
    def setUp(self):
        self.a = RefseqAssembly(
            assembly_accession="GCF_000010525.1",
            ftp_path="ftp://ncbi.gov/abcd",
            asm_name="ASM1052v1",
        )

    def test_parse_summary(self):
        with open(SUMMARY_FP) as f:
            assemblies = list(RefseqAssembly.parse_summary(f))
            self.assertEqual(assemblies[0].accession, "GCF_000010525.1")
            self.assertEqual(
                assemblies[0].ftp_path,
                (
                    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/"
                    "GCF/000/010/525/GCF_000010525.1_ASM1052v1"
                ),
            )

    def test_base_url(self):
        self.assertEqual(self.a.base_url, "https://ncbi.gov/abcd")

    def test_rna_url(self):
        self.assertEqual(
            self.a.rna_url, "https://ncbi.gov/abcd/abcd_rna_from_genomic.fna.gz"
        )

    def test_genome_url(self):
        self.assertEqual(self.a.genome_url, "https://ncbi.gov/abcd/abcd_genomic.fna.gz")
