import pathlib
import shutil
import tempfile
import unittest

from unassign.trimragged import (
    deambiguate, partial_seqs, pairs,
    trim_left, trim_right, main,
    TrimmableSeqs, PrimerMatch,
    PartialMatcher, CompleteMatcher,
)

BSF8 = "AGAGTTTGATCCTGGCTCAG"

class TrimraggedFunctions(unittest.TestCase):
    def test_deambiguate(self):
        self.assertEqual(deambiguate("ACTG"), ["ACTG"])
        self.assertEqual(
            set(deambiguate("CTGCTGCCTYCCGTA")),
            set(["CTGCTGCCTTCCGTA", "CTGCTGCCTCCCGTA"]))

    def test_partial_seqs(self):
        res3 = list(partial_seqs("ABCDEFG", 3))
        self.assertEqual(res3, ["BCDEFG", "CDEFG", "DEFG", "EFG"])

        res4 = list(partial_seqs("ABCDEFG", 4))
        self.assertEqual(res4, ["BCDEFG", "CDEFG", "DEFG"])

class TrimmableSeqsTest(unittest.TestCase):
    def test_get_unmatched_recs(self):
        recs = [
            ("AF403541", "TAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTA"),
            ("AF403542", "GATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGC"),
            ("AF403543", "ACTGTCGTGTCGAAGTGTGGGCGTACGTGTTTGCAACGTGTCAA"),
            ("AF403544", "TAGAGTATGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTA"),
        ]
        s = TrimmableSeqs.from_fasta(MAIN_INPUT.splitlines())
        self.assertEqual(list(s.get_unmatched_recs()), recs)
        class MockMatch(object):
            pass
        s.register_match("AF403544", MockMatch())
        self.assertEqual(list(s.get_unmatched_recs()), recs[0:3])
        self.assertEqual(list(s.get_matched_recs()), recs[3:])

class MatcherFunctions(unittest.TestCase):
    def test_partial_match(self):
        qset = [BSF8]
        seq = "GATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGG"
        trimmed_seq = "GACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGG"
        m = PartialMatcher(qset, 10)
        matchobj = m.find_match(seq)
        self.assertEqual(matchobj.start, 0)
        self.assertEqual(matchobj.end, 13)
        self.assertEqual(trim_left(seq, matchobj), trimmed_seq)

    def test_exact_match(self):
        qset = [BSF8]
        seq = "TAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGG"
        m = CompleteMatcher(qset, max_mismatch=0)
        matchobj = m.find_match(seq)
        self.assertEqual(matchobj.start, 1)
        self.assertEqual(matchobj.end, 21)

    def test_2_mismatches(self):
        qset = [BSF8]
        seq = "TAGAGTAAGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGG"
        m = CompleteMatcher(qset, max_mismatch=2)
        matchobj = m.find_match(seq)
        self.assertEqual(matchobj.start, 1)
        self.assertEqual(matchobj.end, 21)

    def test_trim_right(self):
        seq = "TCCTAGAG"
        class MockMatch(object):
            start = 4
            end = 6
        matchobj = MockMatch()
        self.assertEqual(trim_right(seq, matchobj), "TCCT")

    def test_pairs(self):
        self.assertEqual(
            list(pairs("ABCDEFG")),
            [("A", "B"), ("C", "D"), ("E", "F")])
        
class TrimraggedMain(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_trimragged_main(self):
        input_fp = pathlib.Path(self.test_dir, "in.fasta")
        with input_fp.open("w") as f:
            f.write(MAIN_INPUT)
        output_fp = pathlib.Path(self.test_dir, "out.fasta")
        stats_fp = pathlib.Path(self.test_dir, "stats.txt")
        args = [
            "--input_file", str(input_fp),
            "--trimmed_output_file", str(output_fp),
            "--stats_output_file", str(stats_fp),
            "--query", BSF8,
        ]
        main(args)

        with output_fp.open() as f:
            output_contents = f.read()
        self.assertEqual(output_contents, MAIN_OUTPUT)

        with stats_fp.open() as f:
            output_contents = f.read()
        self.assertEqual(output_contents, MAIN_STATS)

MAIN_INPUT = """\
>AF403541 full BSF8
TAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTA
>AF403542 partial BSF8
GATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGC
>AF403543 nomatch
ACTGTCGTGTCGAAGTGTGGGCGTACGTGTTTGCAACGTGTCAA
>AF403544 BSF8 with mismatch
TAGAGTATGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTA
"""

MAIN_OUTPUT = """\
>AF403541 full BSF8
GACGAACGCTGGCGGCGTGCTTA
>AF403544 BSF8 with mismatch
GACGAACGCTGGCGGCGTGCTTA
>AF403542 partial BSF8
GACGAACGCTGGCGGCGTGCTTAACACATGC
"""

MAIN_STATS = """\
AF403541	Exact	1	21
AF403544	Complete, 1 mismatch	1	21
AF403542	Partial	0	13
AF403543	Unmatched	NA	NA
"""
