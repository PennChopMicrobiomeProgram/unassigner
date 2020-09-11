import gzip
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
            self.assertEqual(
                vals[:5], ["query0", "Aname for2", "2", "0", "265"])
            unassignment_prob = float(vals[5])
            self.assertAlmostEqual(unassignment_prob, 0.0001, 4)

    def test_mismatch_db(self):
        query_fp = os.path.join(self.dir, "query.fa")
        with open(query_fp, "w") as f:
            f.write(query_fasta)

        mismatch_fp = os.path.join(self.dir, "mismatch.txt.gz")
        with gzip.open(mismatch_fp, "wt") as f:
            f.write(mismatch_positions)

        main([
            query_fp,
            "--output_dir", self.dir,
            "--type_strain_fasta", REF_FP,
            "--ref_mismatch_positions", mismatch_fp,
        ])

        with(open(os.path.join(self.dir, "unassigner_output.tsv"))) as f:
            header_line = next(f)
            self.assertTrue(header_line.startswith("query_id"))
            results_line = next(f)
            vals = results_line.rstrip("\n").split("\t")
            self.assertEqual(
                vals[:5], ["query0", "Aname for2", "2", "0", "265"])
            unassignment_prob = float(vals[5])
            # unassignment_prob goes up because we expect to see more
            # mismatches outside the aligned region
            self.assertAlmostEqual(unassignment_prob, 0.0005, 4)

# Exact match to referece sequence 2
query_fasta = """\
>query0
AATGAACGCTGGCGGCGTGCCTAACACATGCAAGTCGTACGAGAAATCCCGAGCTTGCTTGGGAAAGTAAAGTGGCGCACGGGTGAGTAACGCGTGGGTAACCCACCCCCGAATTCGGGATAACTCCGCGAAAGCGGTGCTAATACCGGATAAGACCCCTACCGCTTCGGCGGCAGAGGTAAAAGCTGACCTCTCCATGGAAGTTAGCGTTTGGGGACGGGCCCGCGTCCTATCAGCTTGTTGGTGGGGTAACAGCCCACCAAGG
"""

mismatch_positions = """\
2	ref1	500	501	502	503	504	505
2	ref2	500	501	502	503	504	505
2	ref3	500	501	502	503	504	505
2	ref4	500	501	502	503	504	505
2	ref5	500	501	502	503	504	505
2	ref6	500	501	502	503	504	505
2	ref7	500	501	502	503	504	505
2	ref8	500	501	502	503	504	505
2	ref9	500	501	502	503	504	505
2	ref10	500	501	502	503	504	505
"""
