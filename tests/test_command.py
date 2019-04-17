import os
import unittest
import tempfile

from unassigner.command import main

DATA_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "data")
REF_FP = os.path.join(DATA_DIR, "gg10.fasta")

class CommandTests(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def test_main(self):
        query_fp = os.path.join(self.dir, "query.fa")
        with open(query_fp, "w") as f:
            f.write(query_fasta)
        main([
            query_fp,
            "--output_dir", self.dir,
            "--type_strain_fasta", REF_FP])

        with(open(os.path.join(self.dir, "unassigner_output.tsv"))) as f:
            header_line = next(f)
            self.assertTrue(header_line.startswith("query_id"))
            results_line = next(f)
            vals = results_line.rstrip("\n").split("\t")
            self.assertEqual(vals[0], "query0")
            self.assertEqual(vals[1], "Aname for2")
            self.assertEqual(vals[2], "2")
            unassignment_prob = float(vals[3])
            self.assertLess(unassignment_prob, 0.001)

# Exact match to referece sequence 2
query_fasta = """\
>query0
AATGAACGCTGGCGGCGTGCCTAACACATGCAAGTCGTACGAGAAATCCCGAGCTTGCTTGGGAAAGTAAAGTGGCGCACGGGTGAGTAACGCGTGGGTAACCCACCCCCGAATTCGGGATAACTCCGCGAAAGCGGTGCTAATACCGGATAAGACCCCTACCGCTTCGGCGGCAGAGGTAAAAGCTGACCTCTCCATGGAAGTTAGCGTTTGGGGACGGGCCCGCGTCCTATCAGCTTGTTGGTGGGGTAACAGCCCACCAAGG
"""
