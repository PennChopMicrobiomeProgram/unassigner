import os
import unittest
import tempfile

from unassign.command import main

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
        main([query_fp, self.dir, "-t", REF_FP])

        with(open(os.path.join(self.dir, "unassigner_output.tsv"))) as f:
            header_line = next(f)
            self.assertTrue(header_line.startswith("QueryID"))
            results_line = next(f)
            vals = results_line.rstrip("\n").split("\t")
            match_vals = vals[0:2]
            self.assertEqual(match_vals, ["query0", "2"])
            count_vals = [int(x) for x in vals[2:6]]
            self.assertEqual(count_vals, [0, 265, 1265, 38])
            unassignment_prob = float(vals[6])
            self.assertLess(unassignment_prob, 0.001)

# Exact match to referece sequence 2
query_fasta = """\
>query0
AATGAACGCTGGCGGCGTGCCTAACACATGCAAGTCGTACGAGAAATCCCGAGCTTGCTTGGGAAAGTAAAGTGGCGCACGGGTGAGTAACGCGTGGGTAACCCACCCCCGAATTCGGGATAACTCCGCGAAAGCGGTGCTAATACCGGATAAGACCCCTACCGCTTCGGCGGCAGAGGTAAAAGCTGACCTCTCCATGGAAGTTAGCGTTTGGGGACGGGCCCGCGTCCTATCAGCTTGTTGGTGGGGTAACAGCCCACCAAGG
"""
