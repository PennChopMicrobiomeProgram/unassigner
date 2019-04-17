import pathlib
import shutil
import tempfile
import unittest

from unassigner.trim import (
    deambiguate, partial_seqs,
    trim_left, trim_right, main,
    TrimmableSeqs, PrimerMatch,
    PartialMatcher, CompleteMatcher,
    AlignmentMatcher,
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
    def setUp(self):
        self.recs = [
            ("AF403541", "TAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTA"),
            ("AF403542", "GATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGC"),
            ("AF403543", "ACTGTCGTGTCGAAGTGTGGGCGTACGTGTTTGCAACGTGTCAA"),
            ("AF403544", "TAGAGTATGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTA"),
            ("AF403545", "TCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAACGTGCA"),
        ]
        class MockMatch(object):
            offset = 0
        self.matchobj = MockMatch()

    def test_from_fasta(self):
        s = TrimmableSeqs.from_fasta(MAIN_INPUT.splitlines())
        self.assertEqual(list(s.get_unmatched_recs()), self.recs)

    def test_all_matched(self):
        s = TrimmableSeqs(self.recs)

        # register matches for all but the last item
        for seq_id, _ in self.recs[:-1]:
            s.register_match(seq_id, self.matchobj)
        self.assertFalse(s.all_matched())

        # register a match for the last item
        for seq_id, _ in self.recs[-1:]:
            s.register_match(seq_id, self.matchobj)
        self.assertTrue(s.all_matched())

    def test_get_matched_unmatched_recs(self):
        s = TrimmableSeqs(self.recs)
        s.register_match("AF403541", self.matchobj)
        self.assertEqual(list(s.get_unmatched_recs()), self.recs[1:])
        self.assertEqual(list(s.get_matched_offset0()), self.recs[:1])

    def test_get_replicate_ids_recs(self):
        s = TrimmableSeqs.from_fasta(MAIN_INPUT.splitlines())
        self.assertEqual(
            list(s.get_replicate_ids("AF403541")), ["AF403541", "AF403541b"])
        rep_recs = list(s.get_replicate_recs("AF403541"))
        rep_rec_ids, rep_rec_seqs = zip(*rep_recs)
        self.assertEqual(rep_rec_ids, ("AF403541", "AF403541b"))
        self.assertEqual(rep_rec_seqs[0], rep_rec_seqs[1])

    def test_get_desc(self):
        s = TrimmableSeqs.from_fasta(MAIN_INPUT.splitlines())
        self.assertEqual(s.get_desc("AF403541"), "AF403541 full BSF8")


class CompleteMatcherTests(unittest.TestCase):
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

class PartialMatcherTests(unittest.TestCase):
    def test_partial_match(self):
        qset = [BSF8]
        seq = "GATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGG"
        trimmed_seq =      "GACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGG"
        m = PartialMatcher(qset, 10)
        matchobj = m.find_match(seq)
        self.assertEqual(matchobj.start, 0)
        self.assertEqual(matchobj.end, 13)
        self.assertEqual(trim_left(seq, matchobj), trimmed_seq)


def mock_trimmable_seqs(sseq, qseq, primer_start, primer_end):
    class MockSeqs:
        matches = {
            "A": PrimerMatch(primer_start, primer_end, 0, "Test"),
        }
        def all_matched(self):
            return False
        def get_matched_offset0(self):
            yield ("A", sseq)
        def get_unmatched_recs(self):
            yield ("B", qseq)
    return MockSeqs()


class AlignmentMatcherTests(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_alignment_match_subj_left(self):
        s = mock_trimmable_seqs(
            "TCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGG",
            #|||||||||||
                    "CAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGG",
            0, 11,
        )
        m = AlignmentMatcher(self.test_dir)
        alignment_matches = list(m.find_in_seqs(s))
        match_id, matchobj = alignment_matches[0]
        self.assertEqual(matchobj.start, 0)
        self.assertEqual(matchobj.end, 3)

    def test_alignment_match_middle(self):
        s = mock_trimmable_seqs(
            "TCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGG",
            #          |||||
               "TGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGG",
            10, 15,
        )
        m = AlignmentMatcher(self.test_dir)
        alignment_matches = list(m.find_in_seqs(s))
        match_id, matchobj = alignment_matches[0]
        self.assertEqual((matchobj.start, matchobj.end), (7, 12))

    def test_alignment_match_middle_gaps(self):
        s = mock_trimmable_seqs(
            "TCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGG",
            #          ||//
               "TGGCTCAGGCGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGG",
            10, 15,
        )
        m = AlignmentMatcher(self.test_dir)
        alignment_matches = list(m.find_in_seqs(s))
        match_id, matchobj = alignment_matches[0]
        self.assertEqual((matchobj.start, matchobj.end), (7, 11))

    def test_trim_right(self):
        seq = "TCCTAGAG"
        matchobj = PrimerMatch(4, 6, 0, "")
        self.assertEqual(trim_right(seq, matchobj), "TCCT")


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
            BSF8,
            "--input_file", str(input_fp),
            "--trimmed_output_file", str(output_fp),
            "--stats_output_file", str(stats_fp),
            "--min_partial", "5",
            "--max_mismatch", "1",
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
>AF403541b full BSF8
TAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTA
>AF403542 partial BSF8
GATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGC
>AF403543 nomatch
ACTGTCGTGTCGAAGTGTGGGCGTACGTGTTTGCAACGTGTCAA
>AF403544 BSF8 with mismatch
TAGAGTATGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTA
>AF403545 only 4bp BSF8
TCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAACGTGCA
"""

MAIN_OUTPUT = """\
>AF403541 full BSF8
GACGAACGCTGGCGGCGTGCTTA
>AF403541b full BSF8
GACGAACGCTGGCGGCGTGCTTA
>AF403544 BSF8 with mismatch
GACGAACGCTGGCGGCGTGCTTA
>AF403542 partial BSF8
GACGAACGCTGGCGGCGTGCTTAACACATGC
>AF403545 only 4bp BSF8
GACGAACGCTGGCGGCGTGCTTAACACATGCAAACGTGCA
"""

MAIN_STATS = """\
AF403541	Exact	1	21	0	AGAGTTTGATCCTGGCTCAG
AF403541b	Exact	1	21	0	AGAGTTTGATCCTGGCTCAG
AF403544	Complete, 1 mismatch	1	21	0	AGAGTATGATCCTGGCTCAG
AF403542	Partial	0	13	0	GATCCTGGCTCAG
AF403545	Alignment	0	4	0	TCAG
AF403543	Unmatched	NA	NA	NA	
"""
