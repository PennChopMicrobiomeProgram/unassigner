import unittest

from unassigner.extend_hit import HitExtender

class HitExtenderTests(unittest.TestCase):
    def test_is_global(self):
        hit = {
            "qstart": 1, "qend": 37, "qlen": 37,
            "sstart": 1, "send": 37, "slen": 37,
        }
        self.assertTrue(HitExtender._is_global(hit))

    def test_detect_leftgap(self):
        hit = {
            "qstart": 5, "qend": 28, "qlen": 28,
            "sstart": 45, "send": 68, "slen": 1336,
        }
        self.assertFalse(HitExtender._is_semiglobal(hit))
        hit["qstart"] = 1
        self.assertTrue(HitExtender._is_semiglobal(hit))

    def test_detect_rightgap(self):
        hit = {
            "qstart": 1, "qend": 24, "qlen": 28,
            "sstart": 41, "send": 64, "slen": 1336,
        }
        self.assertFalse(HitExtender._is_semiglobal(hit))
        hit["qend"] = 28
        self.assertTrue(HitExtender._is_semiglobal(hit))

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
