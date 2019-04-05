import os.path
import unittest

from unassign.algorithm import (
    UnassignAligner, BasicAlgorithm,
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

class BasicAlgorithmTests(unittest.TestCase):
    def setUp(self):
        self.ggfp = os.path.join(DATA_DIR, "gg10.fasta")
        a = UnassignAligner(self.ggfp)
        self.basic = BasicAlgorithm(a)

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
