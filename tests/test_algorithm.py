import os.path
import unittest

from unassigner.algorithm import (
    UnassignAligner, UnassignerApp,
    VariableMismatchRate,
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

class VariableMismatchRateTests(unittest.TestCase):
    def setUp(self):
        self.mismatch_fp = os.path.join(DATA_DIR, "mismatch_db.txt")

    def tearDown(self):
        VariableMismatchRate.clear_database()

    def test_load_database(self):
        with open(self.mismatch_fp) as f:
            VariableMismatchRate.load_database(f)

        aab_mismatches = VariableMismatchRate.db["AABF01000111"]
        self.assertEqual(len(aab_mismatches), 3)
        self.assertEqual(aab_mismatches[0], [])
        self.assertEqual(aab_mismatches[1], [173, 234])
        self.assertEqual(aab_mismatches[2], [173, 234, 876, 1268])

        ab0_mismatches = VariableMismatchRate.db["AB004719"]
        self.assertEqual(len(ab0_mismatches), 4)
        self.assertEqual(ab0_mismatches[0], [])
        self.assertEqual(ab0_mismatches[1], [])
        self.assertEqual(ab0_mismatches[2], [
            43, 86, 138,           410, 481, 520, 550, 1388])
        self.assertEqual(ab0_mismatches[3], [
            43,     138, 168, 295, 410, 481, 520, 550, 1388])

        notindb_mismatches = VariableMismatchRate.db["notindb"]
        self.assertEqual(notindb_mismatches, [])

    def test_get_mismatches(self):
        VariableMismatchRate.db["abc"] = [
            [1, 3, 9, 12, 15, 16],
            [2, 8, 10, 11],
            [20],
        ]
        mms = VariableMismatchRate._get_mismatches("abc", 3, 10)
        self.assertEqual(list(mms), [(2, 4), (2, 2), (0, 1)])

class UnassignerAppTests(unittest.TestCase):
    def setUp(self):
        self.ggfp = os.path.join(DATA_DIR, "gg10.fasta")
        a = UnassignAligner(self.ggfp)
        self.app = UnassignerApp(a)

    def test_threshold(self):
        ref_ids = set(str(x) for x in range(1, 10))
        seqs = [
            ("a", "CTTGCTCTCGGGTGACGAGCGGCGGACGGGTGAGTAAT"),
            ("b", "GCGTGGCGAACGGCTGACGAACACGTGG"),
            ]
        all_results = self.app.unassign(seqs)
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
        all_results = self.app.unassign(seqs)
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
        all_results = self.app.unassign(seqs)
        query_id, query_results = next(all_results)
        self.assertEqual(query_id, "query10")
        self.assertEqual(len(query_results), 2)
        top_match = query_results[0]
        self.assertEqual(top_match["typestrain_id"], "10")
        # Expect very low probability of unassignment
        self.assertLess(top_match["probability_incompatible"], 0.001)
