from collections import namedtuple
import os.path
import tempfile
import unittest

from unassign.download import make_blast_db
from unassign.search_blast import (
    UnassignAligner, SemiGlobalAlignment, Alignment, polish_blast_alignment,
    )

DATA_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "data")


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
        a = Alignment.from_blast_hit(self.hit)
        self.assertEqual(a.query_seq, "CCCGGTCCGGTTATT")
        self.assertEqual(a.subject_seq, "CCCGGTCCGGTTAAC")
        self.assertEqual(a.start_idx(a.subject_seq, a.query_seq), 0)
        self.assertEqual(a.end_idx(a.subject_seq, a.query_seq), 15)
        self.assertEqual(a.query_len, 15)
        self.assertEqual(a.subject_len, 15)
        self.assertEqual(a.count_matches(), (13,15))
        
    def test_midgaps(self):
        self.hit.update({
            "qseq":"CCCGGTC--CGGTTATT", "sseq":"CCCGGTCAACGGTTAAC",
            "send":17, "slen":17
            })
        a = Alignment.from_blast_hit(self.hit)
        self.assertEqual(a.query_seq, "CCCGGTC--CGGTTATT")
        self.assertEqual(a.subject_seq, "CCCGGTCAACGGTTAAC")
        self.assertEqual(a.end_idx(a.subject_seq, a.query_seq), 17)
        self.assertEqual(a.query_len, 15)
        self.assertEqual(a.subject_len, 17)
        self.assertEqual(a.count_matches(), (13,17))

class SemiGlobalAlignmentTests(unittest.TestCase):
    def setUp(self):
        self.query_id = "a"
        self.subject_id = "b"
        self.qseq_orj = "CCCGGTCCGGTTATT"
        self.sseq_orj = "CCCGGTCCGGTTAAC"

    def test_no_endgaps(self):
        a = SemiGlobalAlignment(self.query_id, self.qseq_orj, self.subject_id, self.sseq_orj)
        self.assertEqual(a.query_seq,   "CCCGGTCCGGTTATT")
        self.assertEqual(a.subject_seq, "CCCGGTCCGGTTAAC")
        self.assertEqual(a.start_idx(a.subject_seq, a.query_seq), 0)
        self.assertEqual(a.end_idx(a.subject_seq, a.query_seq), 15)
        self.assertEqual(a.query_len, 15)
        self.assertEqual(a.subject_len, 15)
        self.assertEqual(a.count_matches(), (13,15))

    def test_endgaps_both_sides_query(self):
        self.sseq_orj = "CTCAGGATGAACGCTAGCTACAGGCTTAACACATGCAAGT"
        self.qseq_orj =      "GATGAACGCTAGCTTCAGGCTTAAC"
        a = SemiGlobalAlignment(self.query_id, self.qseq_orj, self.subject_id, self.sseq_orj)
        self.assertEqual(a.query_seq, "GATGAACGCTAGCTTCAGGCTTAAC")
        self.assertEqual(a.subject_seq, "GATGAACGCTAGCTACAGGCTTAAC")
        self.assertEqual(a.start_idx(a.subject_seq, a.query_seq), 0)
        self.assertEqual(a.end_idx(a.subject_seq, a.query_seq), 25)
        self.assertEqual(a.query_len, 25)
        self.assertEqual(a.subject_len, 40)
        self.assertEqual(a.count_matches(), (24,25))
        
    def test_endgaps_left_query(self):
        self.sseq_orj = "CTCAGGATGAACGCTAGCTACAGGCTTAAC"
        self.qseq_orj =      "GATGAACGCTAGCTTCAGGCTTAAC"
        a = SemiGlobalAlignment(self.query_id, self.qseq_orj, self.subject_id, self.sseq_orj)
        self.assertEqual(a.query_seq, "GATGAACGCTAGCTTCAGGCTTAAC")
        self.assertEqual(a.subject_seq, "GATGAACGCTAGCTACAGGCTTAAC")
        self.assertEqual(a.start_idx(a.subject_seq, a.query_seq), 0)
        self.assertEqual(a.end_idx(a.subject_seq, a.query_seq), 25)
        self.assertEqual(a.query_len, 25)
        self.assertEqual(a.subject_len, 30)
        self.assertEqual(a.count_matches(), (24,25))

    def test_endgaps_right_query(self):
        self.sseq_orj = "GATGAACGCTAGCTACAGGCTTAACACATGCAAGT"
        self.qseq_orj = "GATGAACGCTAGCTTCAGGCTTAAC"
        a = SemiGlobalAlignment(self.query_id, self.qseq_orj, self.subject_id, self.sseq_orj)
        self.assertEqual(a.query_seq, "GATGAACGCTAGCTTCAGGCTTAAC")
        self.assertEqual(a.subject_seq, "GATGAACGCTAGCTACAGGCTTAAC")
        self.assertEqual(a.start_idx(a.subject_seq, a.query_seq), 0)
        self.assertEqual(a.end_idx(a.subject_seq, a.query_seq), 25)
        self.assertEqual(a.query_len, 25)
        self.assertEqual(a.subject_len, 35)
        self.assertEqual(a.count_matches(), (24,25))

    def test_ragged_left_subject(self):
        self.qseq_orj = "CTCAGGATGAACGCTAGCTACAGGCTTAAC"
        self.sseq_orj =      "GATGAACGCTAGCTTCAGGCTTAACACATG"
        a = SemiGlobalAlignment(self.query_id, self.qseq_orj, self.subject_id, self.sseq_orj)
        self.assertEqual(a.subject_seq, "GATGAACGCTAGCTTCAGGCTTAAC")
        self.assertEqual(a.query_seq, "GATGAACGCTAGCTACAGGCTTAAC")
        self.assertEqual(a.start_idx(a.subject_seq, a.query_seq), 0)
        self.assertEqual(a.end_idx(a.subject_seq, a.query_seq), 25)
        self.assertEqual(a.query_len, 30)
        self.assertEqual(a.subject_len, 30)
        self.assertEqual(a.count_matches(), (24,25))
        
    def test_ragged_right_subject(self):
        self.qseq_orj =      "GATGAACGCTAGCTACAGGCTTAACACATG"
        self.sseq_orj = "CTCAGGATGAACGCTAGCTTCAGGCTTAAC"
        a = SemiGlobalAlignment(self.query_id, self.qseq_orj, self.subject_id, self.sseq_orj)
        self.assertEqual(a.subject_seq, "GATGAACGCTAGCTTCAGGCTTAAC")
        self.assertEqual(a.query_seq, "GATGAACGCTAGCTACAGGCTTAAC")
        self.assertEqual(a.start_idx(a.subject_seq, a.query_seq), 0)
        self.assertEqual(a.end_idx(a.subject_seq, a.query_seq), 25)
        self.assertEqual(a.query_len, 30)
        self.assertEqual(a.subject_len, 30)
        self.assertEqual(a.count_matches(), (24,25))
        
    def test_funky_query_left(self):
        self.qseq_orj =  "TTTTGATGAACGCTAGCTACAGGCTTA"
        self.sseq_orj = "CTCAGGATGAACGCTAGCTTCAGGCTTAAC"
        a = SemiGlobalAlignment(self.query_id, self.qseq_orj, self.subject_id, self.sseq_orj)
        self.assertEqual(a.subject_seq, "TCAGGATGAACGCTAGCTTCAGGCTTA")
        self.assertEqual(a.query_seq,   "TTTTGATGAACGCTAGCTACAGGCTTA")
        self.assertEqual(a.start_idx(a.subject_seq, a.query_seq), 0)
        self.assertEqual(a.end_idx(a.subject_seq, a.query_seq), 27)
        self.assertEqual(a.query_len, 27)
        self.assertEqual(a.subject_len, 30)
        self.assertEqual(a.count_matches(), (23,27))
    

class BlastAlignerTests(unittest.TestCase):
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

    def test_polish_alignment_leftgap(self):
        seqs = [("b", "GCGTGGCGAACGGCTGACGAACACGTGG")]
        hit = {
            "qseqid": "b", "sseqid": "5",
            "qseq": "GCGTGGCGAACGGCTGACGAACACGTGG",
            "sseq": "GCGTGGCGAACGGCTGACGAACACGTGG",
            "qstart": 5, "qend": 28, "qlen": 28,
            "sstart": 45, "send": 68, "slen": 1336,
        }
        self.assertIsInstance(
            polish_blast_alignment(hit, seqs, self.ggfp),
            SemiGlobalAlignment)

    def test_polish_alignment_rightgap(self):
        seqs = [("b", "GCGTGGCGAACGGCTGACGAACACGTGG")]
        hit = {
            "qseqid": "b", "sseqid": "5",
            "qseq": "GCGTGGCGAACGGCTGACGAACACGTGG",
            "sseq": "GCGTGGCGAACGGCTGACGAACACGTGG",
            "qstart": 1, "qend": 24, "qlen": 28,
            "sstart": 41, "send": 64, "slen": 1336,
        }
        self.assertIsInstance(
            polish_blast_alignment(hit, seqs, self.ggfp),
            SemiGlobalAlignment)
        
class HitIdentityTests(unittest.TestCase):
    def test_hit_identity_no_endgaps(self):
        # Hit has 15 positions and 2 mismatches (rightmost columns).
        hit = {
            "qseqid": "a", "sseqid": "b",
            "qseq": "CCCGGTCCGGTTATT",
            "sseq": "CCCGGTCCGGTTAAC",
            "qstart": 1, "qend": 15, "qlen": 15,
            "sstart": 1, "send": 15, "slen": 15,
            }
        a = Alignment.from_blast_hit(hit)
        self.assertEqual(a.count_matches(), (13, 15))

    def test_hit_identity_query_gaps(self):
        # Hit has 14 positions and 2 gaps and 2 mismatches (rightmost columns).
        # Query has 1 nt to right of the alignment.
        # Total matches: 14 - 2 = 12
        # Total positions counted: 14 = query length
        hit = {
            "qseqid": "a", "sseqid": "b",
            "qseq": "CCCGGTCCGGTT--TT",
            "sseq": "CCCGGTCCGGTTAACC",
            "qstart": 1, "qend": 14, "qlen": 14,
            "sstart": 1, "send": 15, "slen": 20,
            }
        a = Alignment.from_blast_hit(hit)
        self.assertEqual(a.count_matches(), (12, 16))

if __name__ == "__main__":
    unittest.main()
