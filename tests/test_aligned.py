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

    def test_region_subject_to_query_no_endgaps(self):
        a = AlignedSubjectQuery(
            ("a", "ABCDEF"),
            ("b", "HIJKLM"))
        # In an alignment with no gaps, the query sequence coordinates should
        # always match the subject sequence coordinates
        self.assertEqual(a.region_subject_to_query(0, 3), (0, 3))
        self.assertEqual(a.region_subject_to_query(1, 5), (1, 5))
        self.assertEqual(a.region_subject_to_query(), (0, 6)) # whole sequence

    def test_region_subject_to_query_with_endgaps(self):
        a = AlignedSubjectQuery(
            ("a", "--ABC-EF---"),
            ("b", "HIJKLMNOPQR"))
        self.assertEqual(a.region_subject_to_query(0, 3), (0, 1)) # A in HIJ
        self.assertEqual(a.region_subject_to_query(1, 6), (0, 3)) # ABC in IJKLM
        self.assertEqual(a.region_subject_to_query(), (0, 5)) # ABCEF in subject

    def test_region_subject_to_query_crazy_alignment(self):
        a = AlignedSubjectQuery(
            ("a", "-A-BC-EF---"),
            ("b", "--HI-JK-LMN"))
        self.assertEqual(a.region_subject_to_query(0, 3), (1, 3)) # BC in HIJ
        self.assertEqual(a.region_subject_to_query(1, 4), (1, 4)) # BCE in IJK
        self.assertEqual(a.region_subject_to_query(), (1, 5)) # BCEF in subject

    def test_region_query_to_subject_no_endgaps(self):
        a = AlignedSubjectQuery(
            ("a", "ABCDEF"),
            ("b", "HIJKLM"))
        # In an alignment with no gaps, the query sequence coordinates should
        # always match the subject sequence coordinates
        self.assertEqual(a.region_query_to_subject(0, 3), (0, 3))
        self.assertEqual(a.region_query_to_subject(1, 5), (1, 5))
        self.assertEqual(a.region_query_to_subject(), (0, 6)) # whole sequence

    def test_region_query_to_subject_with_endgaps(self):
        a = AlignedSubjectQuery(
            ("a", "--ABC-EF---"),
            ("b", "HIJKLMNOPQR"))
        self.assertEqual(a.region_query_to_subject(0, 3), (2, 5)) # JKL in ABC
        self.assertEqual(a.region_query_to_subject(2, 5), (4, 8)) # LMNO in CEF
        self.assertEqual(a.region_query_to_subject(), (2, 8)) # JKLMNO in query

    def test_region_query_to_subject_crazy_alignment(self):
        a = AlignedSubjectQuery(
            ("a", "-A-BC-EF---"),
            ("b", "--HI-JK-LMN"))
        self.assertEqual(a.region_query_to_subject(0, 3), (0, 2)) # HI in ABC
        self.assertEqual(a.region_query_to_subject(1, 4), (1, 4)) # IJK in BCE
        self.assertEqual(a.region_query_to_subject(), (0, 4)) # HIJK in query

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
