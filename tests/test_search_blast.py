from collections import namedtuple
import os.path
import tempfile
import unittest

from unassign.download import make_blast_db
from unassign.search_blast import (
    _hit_identity, BlastAligner, BlastAlignment,
    )

DATA_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "data")

BLAST_NOT_INSTALLED = not BlastAligner.is_installed()

class BlastAlignmentTests(unittest.TestCase):
    def setUp(self):
        # Hit has 15 positions and 2 mismatches (rightmost columns).
        self.hit = {
            "qseqid": "a", "sseqid": "b",
            "qseq": "CCCGGTCCGGTTATT",
            "sseq": "CCCGGTCCGGTTAAC",
            "qstart": 1, "qend": 15, "qlen": 15,
            "sstart": 1, "send": 15, "slen": 15,
            }
        self.pairs = list(zip("CCCGGTCCGGTTATT", "CCCGGTCCGGTTAAC"))

    def test_no_endgaps(self):
        a = BlastAlignment(self.hit)
        self.assertEqual(a.alignment_length(), 15)
        self.assertEqual(a.num_matches(), 13)
        self.assertEqual(a.unaligned_length(), 0)

        self.assertEqual(a.query_seq, "CCCGGTCCGGTTATT")
        self.assertEqual(a.subject_seq, "CCCGGTCCGGTTAAC")
        self.assertEqual(a.start_idx, 0)
        self.assertEqual(a.end_idx, 15)
        self.assertEqual(a.get_pairs(), self.pairs)

    def test_equal_endgaps_left(self):
        # Hit has 15 positions and 2 mismatches (rightmost columns).
        # Both query and subject seqs have 5 unaligned bases to the
        # left of the local alignment region.
        self.hit.update({
            "qstart": 6, "qend": 20, "qlen": 20,
            "sstart": 6, "send": 20, "slen": 20,
            })
        a = BlastAlignment(self.hit)
        self.assertEqual(a.alignment_length(), 20)
        self.assertEqual(a.num_matches(), 13)
        self.assertEqual(a.unaligned_length(), 0)

        self.assertEqual(a.query_seq, "XXXXXCCCGGTCCGGTTATT")
        self.assertEqual(a.subject_seq, "HHHHHCCCGGTCCGGTTAAC")
        self.assertEqual(a.start_idx, 0)
        self.assertEqual(a.end_idx, 20)
        self.assertEqual(a.get_pairs(), list(zip("XXXXX", "HHHHH")) + self.pairs)

    def test_equal_endgaps_right(self):
        # Hit has 15 positions and 2 mismatches (rightmost columns).
        # Both query and subject seqs have 5 unaligned bases to the
        # right of the local alignment region.
        self.hit.update({"qlen": 20, "slen": 20})
        a = BlastAlignment(self.hit)
        self.assertEqual(a.alignment_length(), 20)
        self.assertEqual(a.num_matches(), 13)
        self.assertEqual(a.unaligned_length(), 0)

        self.assertEqual(a.query_seq, "CCCGGTCCGGTTATTXXXXX")
        self.assertEqual(a.subject_seq, "CCCGGTCCGGTTAACHHHHH")
        self.assertEqual(a.start_idx, 0)
        self.assertEqual(a.end_idx, 20)
        self.assertEqual(a.get_pairs(), self.pairs + list(zip("XXXXX", "HHHHH")))

    def test_query_endgaps_left(self):
        # Hit has 15 positions and 2 mismatches (rightmost columns).
        # Query has 3 unaligned bases to the left of the local
        # alignment region.
        self.hit.update({
            "qstart": 4, "qend": 18, "qlen": 18,
            "sstart": 6, "send": 20, "slen": 20,
            })
        a = BlastAlignment(self.hit)
        self.assertEqual(a.alignment_length(), 18)
        self.assertEqual(a.num_matches(), 13)
        self.assertEqual(a.unaligned_length(), 2)

        self.assertEqual(a.query_seq,   "--XXXCCCGGTCCGGTTATT")
        self.assertEqual(a.subject_seq, "HHHHHCCCGGTCCGGTTAAC")
        self.assertEqual(a.start_idx, 2)
        self.assertEqual(a.end_idx, 20)
        self.assertEqual(a.get_pairs(), list(zip("XXX", "HHH")) + self.pairs)

    def test_query_endgaps_right(self):
        # Hit has 15 positions and 2 mismatches (rightmost columns).
        # Query has 3 unaligned bases to the right of the local
        # alignment region.
        self.hit.update({
            "qend": 15, "qlen": 18,
            "send": 15, "slen": 20,
            })
        a = BlastAlignment(self.hit)
        self.assertEqual(a.alignment_length(), 18)
        self.assertEqual(a.num_matches(), 13)
        self.assertEqual(a.unaligned_length(), 2)

        self.assertEqual(a.query_seq,   "CCCGGTCCGGTTATTXXX--")
        self.assertEqual(a.subject_seq, "CCCGGTCCGGTTAACHHHHH")

    def test_subject_endgaps_left(self):
        # Hit has 15 positions and 2 mismatches (rightmost columns).
        # Subject has 3 unaligned bases to the left of the local
        # alignment region.
        self.hit.update({
            "qstart": 6, "qend": 20, "qlen": 20,
            "sstart": 4, "send": 18, "slen": 18,
            })
        a = BlastAlignment(self.hit)
        self.assertEqual(a.alignment_length(), 18)
        self.assertEqual(a.num_matches(), 13)
        self.assertEqual(a.unaligned_length(), 0)

        self.assertEqual(a.query_seq,   "XXXXXCCCGGTCCGGTTATT")
        self.assertEqual(a.subject_seq, "--HHHCCCGGTCCGGTTAAC")

    def test_subject_endgaps_right(self):
        # Hit has 15 positions and 2 mismatches (rightmost columns).
        # Query has 3 unaligned bases to the right of the local
        # alignment region.
        self.hit.update({
            "qstart": 1, "qend": 15, "qlen": 20,
            "sstart": 1, "send": 15, "slen": 18,
            })
        a = BlastAlignment(self.hit)
        self.assertEqual(a.alignment_length(), 18)
        self.assertEqual(a.num_matches(), 13)
        self.assertEqual(a.unaligned_length(), 0)

        self.assertEqual(a.query_seq,   "CCCGGTCCGGTTATTXXXXX")
        self.assertEqual(a.subject_seq, "CCCGGTCCGGTTAACHHH--")

    def test_query_gaps(self):
        # Hit has 15 positions and 3 mismatches (rightmost columns).
        # Query has 1 nt to right of the alignment.
        # Total matches: 15 - 3 = 12
        # Total positions counted: 15 + 1 = 16
        self.hit.update({
            "qseq": "CCCGGTCCGGTT-TT",
            "sseq": "CCCGGTCCGGTTAAC",
            "qstart": 1, "qend": 14, "qlen": 15,
            "sstart": 1, "send": 15, "slen": 20,
            })
        a = BlastAlignment(self.hit)
        self.assertEqual(a.alignment_length(), 16)
        self.assertEqual(a.num_matches(), 12)
        self.assertEqual(a.unaligned_length(), 4)

        self.assertEqual(a.query_seq,   "CCCGGTCCGGTT-TTX----")
        self.assertEqual(a.subject_seq, "CCCGGTCCGGTTAACHHHHH")


class BlastAlignerTests(unittest.TestCase):
    def setUp(self):
        ggfp = os.path.join(DATA_DIR, "gg10.fasta")
        self.a = BlastAligner(ggfp)

    @unittest.skipIf(BLAST_NOT_INSTALLED, "NCBI BLAST is not installed")
    def test_search_species(self):
        seqs = [
            ("a", "CTTGCTCTCGGGTGACGAGCGGCGGACGGGTGAGTAAT"),
            ("b", "GCGTGGCGAACGGCTGACGAACACGTGG"),
            ]
        hits = self.a.search_species(seqs)
        observed = [(hit.query_id, hit.subject_id) for hit in hits]
        expected = [("a", "8"), ("b", "5")]
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
