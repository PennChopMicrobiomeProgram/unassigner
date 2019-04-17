import os.path
import unittest

from unassigner.algorithm import (
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
        self.algo = ThresholdAlgorithm(a)

    def test_threshold(self):
        ref_ids = set(str(x) for x in range(1, 10))
        seqs = [
            ("a", "CTTGCTCTCGGGTGACGAGCGGCGGACGGGTGAGTAAT"),
            ("b", "GCGTGGCGAACGGCTGACGAACACGTGG"),
            ]
        all_results = self.algo.unassign(seqs)
        first_query_id, first_query_results = next(all_results)
        self.assertEqual(first_query_id, "a")
        first_query_match = first_query_results[0]
        self.assertIn(first_query_match["typestrain_id"], ref_ids)

    def test_low_prob(self):
        # Exact match to part of reference sequence 10
        exact_gg10 = (
            "GGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAGCGGCAGCGCGG"
            "GGCAACCTGGCGGCGAGCGGCGAACGGGTGAGTAACGCGTAGGAATCTACCCAGTAG"
            "CGGGGGATAGCCCGGGGAAACTCGGATTAATACCGCATACGCCCTAAGGGGGAAAGC"
            "AGGGGATCTTCGGACCTTGCACTATTGGAAGAGCCTGCGTTGGATTAGCTAGTTGGT"
            "AGGGTAAAGGCCTACCAAGGCGACGATCCATA")
        seqs = [("query0", exact_gg10)]
        all_results = self.algo.unassign(seqs)
        query_id, query_results = next(all_results)
        top_match = query_results[0]
        self.assertEqual(query_id, "query0")
        self.assertEqual(top_match["typestrain_id"], "10")
        # Expect very low probability of unassignment
        self.assertLess(top_match["probability_incompatible"], 0.001)

    def test_multiple_species_hit(self):
        exact_gg10 = (
            "GGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAGCGGCAGCGCGG"
            "GGCAACCTGGCGGCGAGCGGCGAACGGGTGAGTAACGCGTAGGAATCTACCCAGTAG"
            "CGGGGGATAGCCCGGGGAAACTCGGATTAATACCGCATACGCCCTAAGGGGGAAAGC"
            "AGGGGATCTTCGGACCTTGCACTATTGGAAGAGCCTGCGTTGGATTAGCTAGTTGGT"
            "AGGGTAAAGGCCTACCAAGGCGACGATCCATA")
        seqs = [("query10", exact_gg10)]
        all_results = self.algo.unassign(seqs)
        query_id, query_results = next(all_results)
        self.assertEqual(query_id, "query10")
        self.assertEqual(len(query_results), 2)
        top_match = query_results[0]
        self.assertEqual(top_match["typestrain_id"], "10")
        # Expect very low probability of unassignment
        self.assertLess(top_match["probability_incompatible"], 0.001)
