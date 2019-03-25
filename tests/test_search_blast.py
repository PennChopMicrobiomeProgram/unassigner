from collections import namedtuple
import os.path
import tempfile
import unittest

from unassign.download import make_blast_db
from unassign.search_blast import (
    UnassignAligner, BlastRefiner, align_semiglobal,
    )

DATA_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "data")

class AlignSemiglobalTests(unittest.TestCase):
    def setUp(self):
        self.query_id = "a"
        self.subject_id = "b"
        self.qseq_orj = "CCCGGTCCGGTTATT"
        self.sseq_orj = "CCCGGTCCGGTTAAC"

    def test_no_endgaps(self):
        qseq = "CCCGGTCCGGTTATT"
        sseq = "CCCGGTCCGGTTAAC"
        self.assertEqual(align_semiglobal(qseq, sseq), (qseq, sseq))

    def test_endgaps_both_sides_query(self):
        qseq =      "GATGAACGCTAGCTTCAGGCTTAAC"
        sseq = "CTCAGGATGAACGCTAGCTACAGGCTTAACACATGCAAGT"
        aligned_qseq, aligned_sseq = align_semiglobal(qseq, sseq)
        self.assertEqual(aligned_qseq, "-----" + qseq + "----------")
        self.assertEqual(aligned_sseq, sseq)
        
    def test_endgaps_left_query(self):
        qseq =      "GATGAACGCTAGCTTCAGGCTTAAC"
        sseq = "CTCAGGATGAACGCTAGCTACAGGCTTAAC"
        aligned_qseq, aligned_sseq = align_semiglobal(qseq, sseq)
        self.assertEqual(aligned_qseq, "-----" + qseq)
        self.assertEqual(aligned_sseq, sseq)

    def test_endgaps_right_query(self):
        qseq = "GATGAACGCTAGCTTCAGGCTTAAC"
        sseq = "GATGAACGCTAGCTACAGGCTTAACACATGCAAGT"
        aligned_qseq, aligned_sseq = align_semiglobal(qseq, sseq)
        self.assertEqual(aligned_qseq, qseq + "----------")
        self.assertEqual(aligned_sseq, sseq)

    def test_ragged_left_subject(self):
        qseq = "CTCAGGATGAACGCTAGCTACAGGCTTAAC"
        sseq =      "GATGAACGCTAGCTTCAGGCTTAACACATG"
        aligned_qseq, aligned_sseq = align_semiglobal(qseq, sseq)
        self.assertEqual(aligned_qseq, qseq + "-----")
        self.assertEqual(aligned_sseq, "-----" + sseq)
        
    def test_ragged_right_subject(self):
        qseq =      "GATGAACGCTAGCTACAGGCTTAACACATG"
        sseq = "CTCAGGATGAACGCTAGCTTCAGGCTTAAC"
        aligned_qseq, aligned_sseq = align_semiglobal(qseq, sseq)
        self.assertEqual(aligned_qseq, "-----" + qseq)
        self.assertEqual(aligned_sseq, sseq + "-----")
        
    def test_funky_query_left(self):
        qseq =  "TTTTGATGAACGCTAGCTACAGGCTTA"
        sseq = "CTCAGGATGAACGCTAGCTTCAGGCTTAAC"
        aligned_qseq, aligned_sseq = align_semiglobal(qseq, sseq)
        self.assertEqual(aligned_qseq, "-" + qseq + "--")
        self.assertEqual(aligned_sseq, sseq)


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

class BlastRefinerTests(unittest.TestCase):
    def test_polish_alignment_leftgap(self):
        seqs = [("b", "GCGTGGCGAACGGCTGACGAACACGTGG")]
        hit = {
            "qseqid": "b", "sseqid": "5",
            "qseq": "GCGTGGCGAACGGCTGACGAACACGTGG",
            "sseq": "GCGTGGCGAACGGCTGACGAACACGTGG",
            "qstart": 5, "qend": 28, "qlen": 28,
            "sstart": 45, "send": 68, "slen": 1336,
        }
        self.assertTrue(BlastRefiner._needs_realignment(hit))

    def test_polish_alignment_rightgap(self):
        seqs = [("b", "GCGTGGCGAACGGCTGACGAACACGTGG")]
        hit = {
            "qseqid": "b", "sseqid": "5",
            "qseq": "GCGTGGCGAACGGCTGACGAACACGTGG",
            "sseq": "GCGTGGCGAACGGCTGACGAACACGTGG",
            "qstart": 1, "qend": 24, "qlen": 28,
            "sstart": 41, "send": 64, "slen": 1336,
        }
        self.assertTrue(BlastRefiner._needs_realignment(hit))

    def test_is_global(self):
        hit = {
            "qstart": 1, "qend": 37, "qlen": 37,
            "sstart": 1, "send": 37, "slen": 37,
        }
        self.assertTrue(BlastRefiner._is_global(hit))

    def test_add_endgaps_left(self):
        no_endgaps = {"qstart": 1, "sstart": 1}
        self.assertEqual(
            BlastRefiner._add_endgaps_left(no_endgaps, "GG", "CC"), ("", ""))
        query_hang = {"qstart": 5, "sstart": 1}
        self.assertEqual(
            BlastRefiner._add_endgaps_left(query_hang, "ABCDEFGH", "JKLMNOP"),
            ("ABCD", "----"))
        subject_hang = {"qstart": 1, "sstart": 4}
        self.assertEqual(
            BlastRefiner._add_endgaps_left(subject_hang, "ABCDEFGH", "JKLMNOP"),
            ("---", "JKL"))

    def test_add_endgaps_right(self):
        no_endgaps = {"qend": 28, "qlen": 28, "send": 14, "slen": 14}
        self.assertEqual(
            BlastRefiner._add_endgaps_right(no_endgaps, "GG", "CC"), ("", ""))
        query_hang = {"qend": 5, "qlen": 8, "send": 10, "slen": 10}
        self.assertEqual(
            BlastRefiner._add_endgaps_right(query_hang, "ABCDEFG", "JKLMNOPQRS"),
            ("EFG", "---"))
        subject_hang = {"qend": 9, "qlen": 9, "send": 5, "slen": 7}
        self.assertEqual(
            BlastRefiner._add_endgaps_right(subject_hang, "ABCDEFGH", "JKLMNOP"),
            ("--", "OP"))



if __name__ == "__main__":
    unittest.main()
