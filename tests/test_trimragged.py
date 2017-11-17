import pathlib
import shutil
import tempfile
import unittest

import unassign.trimragged

BSF8 = "AGAGTTTGATCCTGGCTCAG"

class TrimraggedFunctions(unittest.TestCase):
    def test_partial_seqs(self):
        res3 = list(unassign.trimragged.partial_seqs("ABCDEFG", 3))
        self.assertEqual(res3, ["BCDEFG", "CDEFG", "DEFG", "EFG"])

        res4 = list(unassign.trimragged.partial_seqs("ABCDEFG", 4))
        self.assertEqual(res4, ["BCDEFG", "CDEFG", "DEFG"])

class MatcherFunctions(unittest.TestCase):
    def test_partial_match(self):
        qset = [BSF8]
        m = unassign.trimragged.PartialMatcher(qset, 10)
        rec = unassign.trimragged.SeqRecord(
            "AF403541",
            "GATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGG")
        res = m.find_match(rec)
        self.assertEqual(res.start, 0)
        self.assertEqual(res.end, 13)
        self.assertEqual(
            res.trim_left().seq, "GACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGG")
        
    def test_exact_match(self):
        qset = [BSF8]
        m = unassign.trimragged.CompleteMatcher(qset, max_mismatch=0)
        rec = unassign.trimragged.SeqRecord(
            "AF403541",
            "TAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGG")
        res = m.find_match(rec)
        self.assertEqual(res.start, 1)
        self.assertEqual(res.end, 21)

    def test_2_mismatches(self):
        qset = [BSF8]
        m = unassign.trimragged.CompleteMatcher(qset, max_mismatch=2)
        rec = unassign.trimragged.SeqRecord(
            "AF403541",
            "TAGAGTAAGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGG")
        res = m.find_match(rec)
        self.assertEqual(res.start, 1)
        self.assertEqual(res.end, 21)

    def test_pairs(self):
        self.assertEqual(
            list(unassign.trimragged.pairs("ABCDEFG")),
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
        unassign.trimragged.main(args)

        # with output_fp.open() as f:
        #     output_contents = f.read()
        # self.assertEqual(output_contents, MAIN_OUTPUT)

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
>AF403542 partial BSF8
GACGAACGCTGGCGGCGTGCTTAACACATGC
>AF403544 BSF8 with mismatch
GACGAACGCTGGCGGCGTGCTTA
"""

MAIN_STATS = """\
AF403541	Exact	1	21
AF403544	Complete, 1 mismatch	1	21
AF403542	Partial	0	13
AF403543	Unmatched	NA	NA
"""
