import unittest

from unassign.alignment import (
    AlignedSubjectQuery, count_matches,
)

class AlignedSubjectQueryTests(unittest.TestCase):
    def test_hit_identity_no_endgaps(self):
        a = AlignedSubjectQuery(
            ("a", "CCCGGTCCGGTTATT"),
            #      |||||||||||||xx
            ("b", "CCCGGTCCGGTTAAC"))
        self.assertEqual(count_matches(a.pairs_query()), (13, 15, 15))

    def test_hit_identity_query_gaps(self):
        a = AlignedSubjectQuery(
            ("a", "CCCGGTCCGGTT--TT-----"),
            #      ||||||||||||  xx
            ("b", "CCCGGTCCGGTTAACCGGGTT"))
        self.assertEqual(count_matches(a.pairs_query()), (12, 14, 16))

    def test_pairs_query_no_endgaps(self):
        a = AlignedSubjectQuery(
            ("a", "ABCDEF"),
            ("b", "HIJKLM"))
        self.assertEqual(
            list(a.pairs_query(0, 3)),
            [("A", "H"), ("B", "I"), ("C", "J")])
        self.assertEqual(
            list(a.pairs_query(1, 5)),
            [("B", "I"), ("C", "J"), ("D", "K"), ("E", "L")])
        self.assertEqual(
            list(a.pairs_query()),
            [
                ("A", "H"), ("B", "I"), ("C", "J"),
                ("D", "K"), ("E", "L"), ("F", "M"),
            ])

    def test_pairs_query_with_endgaps(self):
        a = AlignedSubjectQuery(
            ("a", "--ABC-EF---"),
            ("b", "HIJKLMNOPQR"))
        self.assertEqual(
            list(a.pairs_query(0, 3)),
            [("A", "J"), ("B", "K"), ("C", "L")])
        self.assertEqual(
            list(a.pairs_query(1, 4)),
            [("B", "K"), ("C", "L"), ("-", "M"), ("E", "N")])
        self.assertEqual(
            list(a.pairs_query()),
            [
                ("A", "J"), ("B", "K"), ("C", "L"),
                ("-", "M"), ("E", "N"), ("F", "O"),
            ])

    def test_pairs_query_crazy_alignment(self):
        a = AlignedSubjectQuery(
            ("a", "-A-BC-EF---"),
            ("b", "--HI-JK-LMN"))
        self.assertEqual(
            list(a.pairs_query(0, 3)),
            [("A", "-"), ("-", "H"), ("B", "I"), ("C", "-")])
        self.assertEqual(
            list(a.pairs_query(1, 4)),
            [("B", "I"), ("C", "-"), ("-", "J"), ("E", "K")])
        self.assertEqual(
            list(a.pairs_query()),
            [
                ("A", "-"), ("-", "H"), ("B", "I"), ("C", "-"),
                ("-", "J"), ("E", "K"), ("F", "-"),
            ])

    def test_pairs_subject_no_endgaps(self):
        a = AlignedSubjectQuery(
            ("a", "ABCDEF"),
            ("b", "HIJKLM"))
        self.assertEqual(
            list(a.pairs_subject(0, 3)),
            [("A", "H"), ("B", "I"), ("C", "J")])
        self.assertEqual(
            list(a.pairs_subject(1, 5)),
            [("B", "I"), ("C", "J"), ("D", "K"), ("E", "L")])
        self.assertEqual(
            list(a.pairs_subject()),
            [
                ("A", "H"), ("B", "I"), ("C", "J"),
                ("D", "K"), ("E", "L"), ("F", "M"),
            ])

    def test_pairs_subject_with_endgaps(self):
        a = AlignedSubjectQuery(
            ("a", "--ABC-EF---"),
            ("b", "HIJKLMNOPQR"))
        self.assertEqual(
            list(a.pairs_subject(0, 3)),
            [("-", "H"), ("-", "I"), ("A", "J")])
        self.assertEqual(
            list(a.pairs_subject(1, 4)),
            [("-", "I"), ("A", "J"), ("B", "K")])
        self.assertEqual(
            list(a.pairs_subject()),
            [
                ("-", "H"), ("-", "I"), ("A", "J"),
                ("B", "K"), ("C", "L"), ("-", "M"),
                ("E", "N"), ("F", "O"), ("-", "P"),
                ("-", "Q"), ("-", "R"),
            ])

    def test_pairs_subject_crazy_alignment(self):
        a = AlignedSubjectQuery(
            ("a", "-A-BC-EF---"),
            ("b", "--HI-JK-LMN"))
        self.assertEqual(
            list(a.pairs_subject(0, 3)),
            [("-", "H"), ("B", "I"), ("C", "-"), ("-", "J")])
        self.assertEqual(
            list(a.pairs_subject(1, 4)),
            [("B", "I"), ("C", "-"), ("-", "J"), ("E", "K")])
        self.assertEqual(
            list(a.pairs_subject()),
            [
                ("-", "H"), ("B", "I"), ("C", "-"),
                ("-", "J"), ("E", "K"), ("F", "-"),
                ("-", "L"), ("-", "M"), ("-", "N"),
            ])
