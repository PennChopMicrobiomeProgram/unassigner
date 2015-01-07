import unittest

from unassign.search_blast import hit_identity

class HitIdentityTests(unittest.TestCase):
    def test_hit_identity_no_endgaps(self):
        # Hit has 15 positions and 2 mismatches (rightmost columns).
        hit = {
            "qseq": "CCCGGTCCGGTTATT",
            "sseq": "CCCGGTCCGGTTAAC",
            "qstart": 1, "qend": 15, "qlen": 15,
            "sstart": 1, "send": 15, "slen": 15,
            }
        self.assertEqual(hit_identity(hit), (13, 15))

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
        self.assertEqual(hit_identity(hit), (13, 20))

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
        self.assertEqual(hit_identity(hit), (13, 20))

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
        self.assertEqual(hit_identity(hit), (13, 18))

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
        self.assertEqual(hit_identity(hit), (13, 18))

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
        self.assertEqual(hit_identity(hit), (13, 18))

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
        self.assertEqual(hit_identity(hit), (13, 18))

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
        self.assertEqual(hit_identity(hit), (12, 16))

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
        
        self.assertEqual(hit_identity(hit, 5, 15), (8, 12))


if __name__ == "__main__":
    unittest.main()
