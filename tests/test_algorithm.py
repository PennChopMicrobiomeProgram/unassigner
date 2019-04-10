import os.path
import unittest

from unassign.algorithm import (
    UnassignAligner, ThresholdAlgorithm,
    beta_binomial_pdf, beta_binomial_cdf,
)

DATA_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "data")

class UnassignAlignerTests(unittest.TestCase):
    def setUp(self):
        self.ggfp = os.path.join(DATA_DIR, "gg10.fasta")
        self.a = UnassignAligner(self.ggfp)

    def test_search_species(self):
        seqs = [
            ("a", "CTTGCTCTCGGGTGACGAGCGGCGGACGGGTGAGTAAT"),
            ("b", "GCGTGGCGAACGGCTGACGAACACGTGG"),
            ]
        hits = self.a.search_species(seqs)
        print(hits)
        observed = set((hit.query_id, hit.subject_id) for hit in hits)
        expected = set([("a", "8"), ("b", "5")])
        self.assertEqual(observed, expected)

class FunctionTests(unittest.TestCase):
    def test_beta_binomial(self):
        # These parameters are plotted on the wikipedia page
        self.assertAlmostEqual(beta_binomial_pdf(8, 10, 600, 400), 0.12, 2)
        self.assertAlmostEqual(beta_binomial_cdf(8, 10, 600, 400), 0.95, 2)

    def test_beta_binomial_out_of_bounds(self):
        self.assertEqual(beta_binomial_pdf(15, 10, 10, 10), 0) # k > n
        self.assertEqual(beta_binomial_pdf(3, 0, 10, 10), 0) # n == 0
        self.assertEqual(beta_binomial_pdf(-3, 5, 10, 10), 0) # k < 0


class ThresholdAlgorithmTests(unittest.TestCase):
    def setUp(self):
        self.ggfp = os.path.join(DATA_DIR, "gg10.fasta")
        a = UnassignAligner(self.ggfp)
        self.basic = ThresholdAlgorithm(a)

    def test_basic(self):
        seqs = [
            ("a", "CTTGCTCTCGGGTGACGAGCGGCGGACGGGTGAGTAAT"),
            ("b", "GCGTGGCGAACGGCTGACGAACACGTGG"),
            ]
        res = list(self.basic.unassign(seqs))
        res_query_ids = set(x[0] for x in res)
        self.assertEqual(res_query_ids, set(("a", "b")))
        res_species_ids = set(x[1] for x in res)
        ref_ids = set(str(x) for x in range(1, 10))
        self.assertLessEqual(res_species_ids, ref_ids)

    def test_low_prob(self):
        # Exact match to part of reference sequence 10
        exact_gg10 = (
            "GGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAGCGGCAGCGCGG"
            "GGCAACCTGGCGGCGAGCGGCGAACGGGTGAGTAACGCGTAGGAATCTACCCAGTAG"
            "CGGGGGATAGCCCGGGGAAACTCGGATTAATACCGCATACGCCCTAAGGGGGAAAGC"
            "AGGGGATCTTCGGACCTTGCACTATTGGAAGAGCCTGCGTTGGATTAGCTAGTTGGT"
            "AGGGTAAAGGCCTACCAAGGCGACGATCCATA")
        seqs = [("query0", exact_gg10)]
        res = list(self.basic.unassign(seqs))
        self.assertEqual(res[0][0], "query0")
        self.assertEqual(res[0][1], "10")
        # Expect very low probability of unassignment
        self.assertLess(res[0][6], 0.001)
