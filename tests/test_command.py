import gzip
import os
import unittest
import tempfile

from unassigner.command import main
from unassigner.algorithm import VariableMismatchRate

DATA_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "data")
SPECIES_FP = os.path.join(DATA_DIR, "species.fasta")


class CommandTests(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        VariableMismatchRate.clear_database()

    def test_r_gnavus(self):
        query_fp = os.path.join(self.dir, "query.fa")
        with open(query_fp, "w") as f:
            f.write(closeto_r_gnavus_fasta)
        main([
            query_fp,
            "--output_dir", self.dir,
            "--type_strain_fasta", SPECIES_FP])

        with(open(os.path.join(self.dir, "unassigner_output.tsv"))) as f:
            header_line = next(f)
            self.assertTrue(header_line.startswith("query_id"))
            results_line = next(f)
            vals = results_line.rstrip("\n").split("\t")
            self.assertEqual(vals[0], "closeto_r_gnavus")
            self.assertEqual(vals[1], "Ruminococcus gnavus")
            self.assertEqual(vals[2], "X94967")
            self.assertEqual(vals[3], "2")
            self.assertEqual(vals[4], "320")
            unassignment_prob = float(vals[5])
            self.assertAlmostEqual(unassignment_prob, 0.003, 3)

    def test_mismatch_db(self):
        query_fp = os.path.join(self.dir, "query.fa")
        with open(query_fp, "w") as f:
            f.write(closeto_r_gnavus_fasta)

        mismatch_fp = os.path.join(self.dir, "mismatch.txt.gz")
        with gzip.open(mismatch_fp, "wt") as f:
            f.write(fake_mismatch_positions)

        main([
            query_fp,
            "--output_dir", self.dir,
            "--type_strain_fasta", SPECIES_FP,
            "--ref_mismatch_positions", mismatch_fp,
        ])

        with(open(os.path.join(self.dir, "unassigner_output.tsv"))) as f:
            header_line = next(f)
            self.assertTrue(header_line.startswith("query_id"))
            results_line = next(f)
            vals = results_line.rstrip("\n").split("\t")
            self.assertEqual(vals[0], "closeto_r_gnavus")
            self.assertEqual(vals[1], "Ruminococcus gnavus")
            self.assertEqual(vals[2], "X94967")
            self.assertEqual(vals[3], "2")
            self.assertEqual(vals[4], "320")
            unassignment_prob = float(vals[5])
            self.assertAlmostEqual(unassignment_prob, 0.085, 3)

    def test_threshold(self):
        query_fp = os.path.join(self.dir, "query.fa")
        with open(query_fp, "w") as f:
            f.write(closeto_r_gnavus_fasta)

        main([
            query_fp,
            "--output_dir", self.dir,
            "--threshold", "0.988",
            "--type_strain_fasta", SPECIES_FP,
        ])

        with(open(os.path.join(self.dir, "unassigner_output.tsv"))) as f:
            header_line = next(f)
            self.assertTrue(header_line.startswith("query_id"))
            results_line = next(f)
            vals = results_line.rstrip("\n").split("\t")
            self.assertEqual(vals[0], "closeto_r_gnavus")
            self.assertEqual(vals[1], "Ruminococcus gnavus")
            self.assertEqual(vals[2], "X94967")
            self.assertEqual(vals[3], "2")
            self.assertEqual(vals[4], "320")
            unassignment_prob = float(vals[5])
            self.assertAlmostEqual(unassignment_prob, 0.147, 3)

    def test_soft_threshold(self):
        query_fp = os.path.join(self.dir, "query.fa")
        with open(query_fp, "w") as f:
            f.write(closeto_r_gnavus_fasta)

        main([
            query_fp,
            "--output_dir", self.dir,
            "--soft_threshold",
            "--type_strain_fasta", SPECIES_FP,
        ])

        with(open(os.path.join(self.dir, "unassigner_output.tsv"))) as f:
            header_line = next(f)
            self.assertTrue(header_line.startswith("query_id"))
            results_line = next(f)
            vals = results_line.rstrip("\n").split("\t")
            self.assertEqual(vals[0], "closeto_r_gnavus")
            self.assertEqual(vals[1], "Ruminococcus gnavus")
            self.assertEqual(vals[2], "X94967")
            self.assertEqual(vals[3], "2")
            self.assertEqual(vals[4], "320")
            unassignment_prob = float(vals[5])
            self.assertAlmostEqual(unassignment_prob, 0.4085, 3)


# Reference mismatches only occur outside the aligned region
fake_mismatch_positions = """\
X94967	ref0	500	501	502	503	504	505
X94967	ref1	500	501	502	503	504	505
X94967	ref2	500	501	502	503	504	505
X94967	ref3	500	501	502	503	504	505
X94967	ref4	500	501	502	503	504	505
X94967	ref5	500	501	502	503	504	505
X94967	ref6	500	501	502	503	504	505
X94967	ref7	500	501	502	503	504	505
X94967	ref8	500	501	502	503	504	505
X94967	ref9	500	501	502	503	504	505
"""

closeto_r_gnavus_fasta = """\
>closeto_r_gnavus
TGAACGCTGGCGGCGTGCTTATCACATGCAAGTCGAGCGAAGCACCTTGACGGATTTCTTCGGATTGAAGCCTTGGTGACTGAGCGGCGGACGGGTGAGTAACGCGTGGGTAACCTGCCACATACAGGGGGATAACAGTTGGAAACGNCTGCTAATACCGCATAAGCGCACAGTACCGCATGGTACGGTGTGAAAAACTCCGGTGGTATGAGATGGACCCGCGTCTGATTAGGTAGTTGGTGGGGTAACGGCCTACCAAGCCGACGATCAGTAGCCGACCTGAGAGGGTGACCGGCCACATTGGGACTGAGACACGGCCC
"""
