import unittest

from unassign.alignment import (
    AlignedPair, AlignedRegion, count_matches,
)

class AlignedPairTests(unittest.TestCase):
    def test_hit_identity_no_endgaps(self):
        a = AlignedPair(
            ("a", "CCCGGTCCGGTTATT"),
            #      |||||||||||||xx
            ("b", "CCCGGTCCGGTTAAC"))
        r = AlignedRegion.from_query(a)
        self.assertEqual(count_matches(r.pairs()), (13, 15, 15))

    def test_hit_identity_query_gaps(self):
        a = AlignedPair(
            ("a", "CCCGGTCCGGTT--TT-----"),
            #      ||||||||||||  xx
            ("b", "CCCGGTCCGGTTAACCGGGTT"))
        r = AlignedRegion.from_query(a)
        self.assertEqual(count_matches(r.pairs()), (12, 14, 16))

    def test_pairs_query_no_endgaps(self):
        a = AlignedPair(
            ("a", "ABCDEF"),
            ("b", "HIJKLM"))
        r = AlignedRegion.from_query(a, 0, 3)
        self.assertEqual(
            list(r.pairs()),
            [("A", "H"), ("B", "I"), ("C", "J")])
        r = AlignedRegion.from_query(a, 1, 5)
        self.assertEqual(
            list(r.pairs()),
            [("B", "I"), ("C", "J"), ("D", "K"), ("E", "L")])
        r = AlignedRegion.from_query(a)
        self.assertEqual(
            list(r.pairs()),
            [
                ("A", "H"), ("B", "I"), ("C", "J"),
                ("D", "K"), ("E", "L"), ("F", "M"),
            ])

    def test_pairs_query_with_endgaps(self):
        a = AlignedPair(
            ("a", "--ABC-EF---"),
            ("b", "HIJKLMNOPQR"))
        r = AlignedRegion.from_query(a, 0, 3)
        self.assertEqual(
            list(r.pairs()),
            [("A", "J"), ("B", "K"), ("C", "L")])
        r = AlignedRegion.from_query(a, 1, 4)
        self.assertEqual(
            list(r.pairs()),
            [("B", "K"), ("C", "L"), ("-", "M"), ("E", "N")])
        r = AlignedRegion.from_query(a)
        self.assertEqual(
            list(r.pairs()),
            [
                ("A", "J"), ("B", "K"), ("C", "L"),
                ("-", "M"), ("E", "N"), ("F", "O"),
            ])

    def test_pairs_query_crazy_alignment(self):
        a = AlignedPair(
            ("a", "-A-BC-EF---"),
            ("b", "--HI-JK-LMN"))
        r = AlignedRegion.from_query(a, 0, 3)
        self.assertEqual(
            list(r.pairs()),
            [("A", "-"), ("-", "H"), ("B", "I"), ("C", "-")])
        r = AlignedRegion.from_query(a, 1, 4)
        self.assertEqual(
            list(r.pairs()),
            [("B", "I"), ("C", "-"), ("-", "J"), ("E", "K")])
        r = AlignedRegion.from_query(a)
        self.assertEqual(
            list(r.pairs()),
            [
                ("A", "-"), ("-", "H"), ("B", "I"), ("C", "-"),
                ("-", "J"), ("E", "K"), ("F", "-"),
            ])

    def test_region_subject_to_query_no_endgaps(self):
        a = AlignedPair(
            ("a", "ABCDEF"),
            ("b", "HIJKLM"))
        # In an alignment with no gaps, the query sequence coordinates should
        # always match the subject sequence coordinates
        r = AlignedRegion.from_subject(a, 0, 3)
        self.assertEqual(r.in_alignment(), (0, 3))
        rq = AlignedRegion.from_query(a, 0, 3)
        self.assertEqual(r, rq)

        r = AlignedRegion.from_subject(a, 1, 5)
        self.assertEqual(r.in_alignment(), (1, 5))
        rq = AlignedRegion.from_query(a, 1, 5)
        self.assertEqual(r, rq)

        r = AlignedRegion.from_subject(a)
        self.assertEqual(r.in_alignment(), (0, 6))
        rq = AlignedRegion.from_query(a)
        self.assertEqual(r, rq)

    def test_region_subject_to_query_with_endgaps(self):
        a = AlignedPair(
            ("a", "--ABC-EF---"),
            ("b", "HIJKLMNOPQR"))

        r = AlignedRegion.from_subject(a, 0, 3)
        self.assertEqual(r.in_subject(), (0, 3)) # HIJ
        self.assertEqual(r.in_alignment(), (0, 3))
        self.assertEqual(r.in_query(), (0, 1)) # A in HIJ

        r = AlignedRegion.from_subject(a, 1, 6)
        self.assertEqual(r.in_subject(), (1, 6)) # IJKLM
        self.assertEqual(r.in_alignment(), (1, 6))
        self.assertEqual(r.in_query(), (0, 3)) # ABC in IJKLM

        r = AlignedRegion.from_subject(a)
        self.assertEqual(r.in_subject(), (0, 11)) # whole sequence
        self.assertEqual(r.in_alignment(), (0, 11))
        self.assertEqual(r.in_query(), (0, 5)) # ABCEF in subject

    def test_region_subject_to_query_crazy_alignment(self):
        a = AlignedPair(
            ("a", "-A-BC-EF---"),
            ("b", "--HI-JK-LMN"))

        r = AlignedRegion.from_subject(a, 0, 3)
        self.assertEqual(r.in_subject(), (0, 3)) # HIJ
        self.assertEqual(r.in_alignment(), (2, 6)) # HI-J
        self.assertEqual(r.in_query(), (1, 3)) # BC in HIJ

        r = AlignedRegion.from_subject(a, 1, 4)
        self.assertEqual(r.in_subject(), (1, 4)) # IJK
        self.assertEqual(r.in_alignment(), (3, 7)) # I-JK
        self.assertEqual(r.in_query(), (1, 4)) # BCE in IJK

        r = AlignedRegion.from_subject(a)
        self.assertEqual(r.in_subject(), (0, 7)) # whole sequence, HIJKLMN
        self.assertEqual(r.in_alignment(), (2, 11)) # HI-JK-LMN
        self.assertEqual(r.in_query(), (1, 5)) # BCEF in subject

    def test_region_subject_to_query_00(self):
        a = AlignedPair(
            ("a", "ABCDEFG"),
            ("b", "---KLMN"),
        )

        r = AlignedRegion.from_subject(a, 0, 2)
        self.assertEqual(r.in_subject(), (0, 2)) # KL
        self.assertEqual(r.in_alignment(), (3, 5)) # KL
        self.assertEqual(r.in_query(), (3, 5)) # DE

        r = AlignedRegion.from_subject(a, 0, 0)
        self.assertEqual(r.in_subject(), (0, 0)) # empty sequence
        self.assertEqual(r.in_alignment(), (3, 3)) # --- | KLMN
        self.assertEqual(r.in_query(), (3, 3)) # ABC | DEFG


class AlignedRegionTests(unittest.TestCase):
    def test_from_subject_no_endgaps(self):
        a = AlignedPair(
            ("a", "ABCDEF"),
            ("b", "HIJKLM"))
        r = AlignedRegion.from_subject(a, 2, 5)
        self.assertEqual(r.start_idx, 2)
        self.assertEqual(r.end_idx, 5)

    def test_from_subject_with_endgaps(self):
        a = AlignedPair(
            ("a", "--ABC-EF---"),
            ("b", "HIJKLMNOPQR"))
        r = AlignedRegion.from_subject(a, 1, 7)
        self.assertEqual(r.start_idx, 1)
        self.assertEqual(r.end_idx, 7)

    def test_from_subject_region_crazy(self):
        a = AlignedPair(
            ("a", "-A-BC-EF---"),
            ("b", "--HI-JK-LMN"))
        r = AlignedRegion.from_subject(a, 0, 3)
        self.assertEqual(r.start_idx, 2)
        self.assertEqual(r.end_idx, 6)

    def test_from_query_no_endgaps(self):
        a = AlignedPair(
            ("a", "ABCDEF"),
            ("b", "HIJKLM"))
        r = AlignedRegion.from_query(a, 2, 5)
        self.assertEqual(r.start_idx, 2)
        self.assertEqual(r.end_idx, 5)
        self.assertEqual(r.in_query(), (2, 5))
        self.assertEqual(r.in_subject(), (2, 5))
        self.assertEqual(r.query_offset(), 0)
        self.assertEqual(r.subject_offset(), 0)

    def test_from_query_with_endgaps(self):
        a = AlignedPair(
            ("a", "--ABC-EF---"),
            ("b", "HIJKLMNOPQR"))
        r = AlignedRegion.from_query(a, 1, 4)
        self.assertEqual(r.start_idx, 3)
        self.assertEqual(r.end_idx, 7)
        self.assertEqual(r.in_query(), (1, 4))
        self.assertEqual(r.in_subject(), (3, 7))
        self.assertEqual(r.query_offset(), 0)
        self.assertEqual(r.subject_offset(), 0)

    def test_from_query_region_crazy(self):
        a = AlignedPair(
            ("a", "-A-BC-EF---"),
            #       ||||          (1, 5)
            ("b", "--HI-JK-LMN"))
        r = AlignedRegion.from_query(a, 0, 3)
        self.assertEqual(r.start_idx, 1)
        self.assertEqual(r.end_idx, 5)
        self.assertEqual(r.in_query(), (0, 3))
        self.assertEqual(r.in_subject(), (0, 2))
        self.assertEqual(r.query_offset(), 0)
        self.assertEqual(r.subject_offset(), 0)

    def test_query_offset_right(self):
        a = AlignedPair(
            ("a", "------ABCDEF"),
            #       |||
            ("b", "GHIJKLMNOP--"))
        r = AlignedRegion(a, 1, 4)
        self.assertEqual(r.in_query(), (0, 0))
        self.assertEqual(r.in_subject(), (1, 4))
        self.assertEqual(r.query_offset(), 2)
        self.assertEqual(r.subject_offset(), 0)

    def test_query_offset_left(self):
        a = AlignedPair(
            ("a", "ABCDEF------"),
            #             |||
            ("b", "--JKLMNOPQRS"))
        r = AlignedRegion(a, 7, 10)
        self.assertEqual(r.in_query(), (6, 6))
        self.assertEqual(r.in_subject(), (5, 8))
        self.assertEqual(r.query_offset(), -1)
        self.assertEqual(r.subject_offset(), 0)

    def test_aligned_start_idx(self):
        seq = "---AB-C-DEF--"
        self.assertEqual(AlignedRegion.aligned_start_idx(seq, 0), 3)
        self.assertEqual(AlignedRegion.aligned_start_idx(seq, 1), 4)
        self.assertEqual(AlignedRegion.aligned_start_idx(seq, 2), 6)
        self.assertEqual(AlignedRegion.aligned_start_idx(seq, 6), 11)

    def test_aligned_end_idx(self):
        seq = "---AB-C-DEF--"
        self.assertEqual(AlignedRegion.aligned_end_idx(seq, 0), 3)
        self.assertEqual(AlignedRegion.aligned_end_idx(seq, 1), 4)
        self.assertEqual(AlignedRegion.aligned_end_idx(seq, 2), 5)
        self.assertEqual(AlignedRegion.aligned_start_idx(seq, 6), 11)
