from collections import namedtuple
import os.path
import tempfile
import unittest

from unassign.download import make_blast_db
from unassign.search_blast import _hit_identity, BlastAligner

DATA_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "data")


MockAlignment = namedtuple("MockAlignment", "subject_id")


class BlastAlignerTests(unittest.TestCase):
    def setUp(self):
        ggfp = os.path.join(DATA_DIR, "gg10.fasta")
        if not os.path.exists(ggfp + ".nin"):
            make_blast_db(ggfp)
        self.a = BlastAligner(ggfp, ggfp)

    def test_search_species(self):
        seqs = [
            ("a", "CTTGCTCTCGGGTGACGAGCGGCGGACGGGTGAGTAAT"),
            ("b", "GCGTGGCGAACGGCTGACGAACACGTGG"),
            ]
        hits = self.a.search_species(seqs)
        observed = [(hit.query_id, hit.subject_id) for hit in hits]
        expected = [("a", "8"), ("b", "5")]
        self.assertEqual(observed, expected)

    def test_search_refseqs(self):
        hits = [MockAlignment("1"), MockAlignment("2")]
        hits = self.a.search_refseqs(hits)
        observed = [(hit.query_id, hit.subject_id) for hit in hits]
        expected = [
            ('1', '1'), ('1', '6'), ('1', '7'), ('1', '10'), ('1', '10'),
            ('1', '3'), ('1', '2'), ('1', '8'), ('1', '8'), ('1', '4'),
            ('1', '9'), ('1', '5'), ('2', '2'), ('2', '10'), ('2', '3'),
            ('2', '8'), ('2', '6'), ('2', '1'), ('2', '7'), ('2', '4'),
            ('2', '5'), ('2', '9'), ('2', '9')]
        self.assertEqual(observed, expected)


class HitIdentityTests(unittest.TestCase):
    def test_hit_identity_no_endgaps(self):
        # Hit has 15 positions and 2 mismatches (rightmost columns).
        hit = {
            "qseq": "CCCGGTCCGGTTATT",
            "sseq": "CCCGGTCCGGTTAAC",
            "qstart": 1, "qend": 15, "qlen": 15,
            "sstart": 1, "send": 15, "slen": 15,
            }
        self.assertEqual(_hit_identity(hit), (13, 15))

    def test_hit_identity_equal_endgaps_left(self):
        # Hit has 15 positions and 2 mismatches (rightmost columns).
        # Both query and subject seqs have 5 unaligned bases to the
        # left of the local alignment region.
        hit = {
            "qseq": "CCCGGTCCGGTTATT",
            "sseq": "CCCGGTCCGGTTAAC",
            "qstart": 6, "qend": 20, "qlen": 20,
            "sstart": 6, "send": 20, "slen": 20,
            }
        self.assertEqual(_hit_identity(hit), (13, 20))

    def test_hit_identity_equal_endgaps_right(self):
        # Hit has 15 positions and 2 mismatches (rightmost columns).
        # Both query and subject seqs have 5 unaligned bases to the
        # right of the local alignment region.
        hit = {
            "qseq": "CCCGGTCCGGTTATT",
            "sseq": "CCCGGTCCGGTTAAC",
            "qstart": 1, "qend": 15, "qlen": 20,
            "sstart": 1, "send": 15, "slen": 20,
            }
        self.assertEqual(_hit_identity(hit), (13, 20))

    def test_hit_identity_query_endgaps_left(self):
        # Hit has 15 positions and 2 mismatches (rightmost columns).
        # Query has 3 unaligned bases to the left of the local
        # alignment region.
        hit = {
            "qseq": "CCCGGTCCGGTTATT",
            "sseq": "CCCGGTCCGGTTAAC",
            "qstart": 4, "qend": 18, "qlen": 18,
            "sstart": 6, "send": 20, "slen": 20,
            }
        self.assertEqual(_hit_identity(hit), (13, 18))

    def test_hit_identity_query_endgaps_right(self):
        # Hit has 15 positions and 2 mismatches (rightmost columns).
        # Query has 3 unaligned bases to the right of the local
        # alignment region.
        hit = {
            "qseq": "CCCGGTCCGGTTATT",
            "sseq": "CCCGGTCCGGTTAAC",
            "qstart": 1, "qend": 15, "qlen": 18,
            "sstart": 1, "send": 15, "slen": 20,
            }
        self.assertEqual(_hit_identity(hit), (13, 18))

    def test_hit_identity_subject_endgaps_left(self):
        # Hit has 15 positions and 2 mismatches (rightmost columns).
        # Subject has 3 unaligned bases to the left of the local
        # alignment region.
        hit = {
            "qseq": "CCCGGTCCGGTTATT",
            "sseq": "CCCGGTCCGGTTAAC",
            "qstart": 6, "qend": 20, "qlen": 20,
            "sstart": 4, "send": 18, "slen": 18,
            }
        self.assertEqual(_hit_identity(hit), (13, 18))

    def test_hit_identity_subject_endgaps_right(self):
        # Hit has 15 positions and 2 mismatches (rightmost columns).
        # Query has 3 unaligned bases to the right of the local
        # alignment region.
        hit = {
            "qseq": "CCCGGTCCGGTTATT",
            "sseq": "CCCGGTCCGGTTAAC",
            "qstart": 1, "qend": 15, "qlen": 20,
            "sstart": 1, "send": 15, "slen": 18,
            }
        self.assertEqual(_hit_identity(hit), (13, 18))

    def test_hit_identity_query_gaps(self):
        # Hit has 15 positions and 3 mismatches (rightmost columns).
        # Query has 1 nt to right of the alignment.
        # Total matches: 15 - 3 = 12
        # Total positions counted: 15 + 1 = 16
        hit = {
            "qseq": "CCCGGTCCGGTT-TT",
            "sseq": "CCCGGTCCGGTTAAC",
            "qstart": 1, "qend": 14, "qlen": 15,
            "sstart": 1, "send": 15, "slen": 20,
            }
        self.assertEqual(_hit_identity(hit), (12, 16))

    def test_hit_identity_region(self):
        # Hit has 15 positions and 3 mismatches (rightmost columns).
        # Query has 1 nt to right of the alignment.
        # Specified region begins in position 5 of query (first 4
        #   positions of query ignored).
        # Total matches: 15 - 4 - 3 = 8
        # Total positions counted: 15 - 4 + 1 = 12
        hit = {
            "qseq": "CCCGGTCCGGTT-TT",
            "sseq": "CCCGGTCCGGTTAAC",
            "qstart": 1, "qend": 14, "qlen": 15,
            "sstart": 1, "send": 15, "slen": 20,
            }
        
        self.assertEqual(_hit_identity(hit, 5, 15), (8, 12))


if __name__ == "__main__":
    unittest.main()
