import unittest

from unassigner.align import (
    HitExtender, align_semiglobal,
    )

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


class HitExtenderTests(unittest.TestCase):
    def test_is_global(self):
        hit = {
            "qstart": 1, "qend": 37, "qlen": 37,
            "sstart": 1, "send": 37, "slen": 37,
        }
        self.assertTrue(HitExtender._is_global(hit))

    def test_realign_leftgap(self):
        hit = {
            "qstart": 5, "qend": 28, "qlen": 28,
            "sstart": 45, "send": 68, "slen": 1336,
        }
        self.assertTrue(HitExtender._needs_realignment(hit))
        hit["qstart"] = 1
        self.assertFalse(HitExtender._needs_realignment(hit))

    def test_realign_rightgap(self):
        hit = {
            "qstart": 1, "qend": 24, "qlen": 28,
            "sstart": 41, "send": 64, "slen": 1336,
        }
        self.assertTrue(HitExtender._needs_realignment(hit))
        hit["qend"] = 28
        self.assertFalse(HitExtender._needs_realignment(hit))

    def test_add_endgaps_left(self):
        no_endgaps = {"qstart": 1, "sstart": 1}
        self.assertEqual(
            HitExtender._add_endgaps_left(no_endgaps, "GG", "CC"), ("", ""))
        query_hang = {"qstart": 5, "sstart": 1}
        self.assertEqual(
            HitExtender._add_endgaps_left(query_hang, "ABCDEFGH", "JKLMNOP"),
            ("ABCD", "----"))
        subject_hang = {"qstart": 1, "sstart": 4}
        self.assertEqual(
            HitExtender._add_endgaps_left(subject_hang, "ABCDEFGH", "JKLMNOP"),
            ("---", "JKL"))

    def test_add_endgaps_right(self):
        no_endgaps = {"qend": 28, "qlen": 28, "send": 14, "slen": 14}
        self.assertEqual(
            HitExtender._add_endgaps_right(no_endgaps, "GG", "CC"), ("", ""))
        query_hang = {"qend": 5, "qlen": 8, "send": 10, "slen": 10}
        self.assertEqual(
            HitExtender._add_endgaps_right(query_hang, "ABCDEFG", "JKLMNOPQRS"),
            ("EFG", "---"))
        subject_hang = {"qend": 9, "qlen": 9, "send": 5, "slen": 7}
        self.assertEqual(
            HitExtender._add_endgaps_right(subject_hang, "ABCDEFGH", "JKLMNOP"),
            ("--", "OP"))


if __name__ == "__main__":
    unittest.main()
