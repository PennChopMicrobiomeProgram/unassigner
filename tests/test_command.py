from collections import namedtuple
import unittest

from unassign.algorithm import RefseqsAlgorithm

class MockAlignment(object):
    def __init__(self, query, subj, start, end):
        self.query_id= query
        self.subject_id = subj
        self.start_idx = start
        self.end_idx = end

    def count_matches(self, start=None, end=None):
        return (1, 1)



class NullAligner(object):
    def search_species(self, seqs):
        return []

    def search_refseqs(self, hits):
        return []


class NoRefAligner(NullAligner):
    def search_species(self, seqs):
        return [
            MockAlignment(seq_id, "S1", 1, len(seq))
            for seq_id, seq in seqs]


class SingleRefAligner(NoRefAligner):
    def search_refseqs(self, hits):
        return [
            MockAlignment(hit.subject_id, "R1", 1, 500)
            for hit in hits]


class UnassignerTests(unittest.TestCase):
    def test_unassign(self):
        seqs = [
            ("a", "CTTGCTCTCGGGTGACGAGCGGCGGACGGGTGAGTAAT")
            ]
        u = RefseqsAlgorithm(SingleRefAligner())
        self.assertEqual(
            list(u.unassign(seqs)), [("a", "S1", 1, 1, "R1", 1, 1)])


if __name__ == "__main__":
    unittest.main()
