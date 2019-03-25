import os.path
import unittest

from unassign.algorithm import UnassignAligner

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
        observed = [(hit.query_id, hit.subject_id) for hit in hits]
        expected = [("a", "8"), ("b", "5")]
        self.assertEqual(observed, expected)

