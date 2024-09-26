import tempfile
import unittest

from unassigner.parse import parse_fasta, parse_results, load_fasta, write_fasta


class FastaTests(unittest.TestCase):
    def test_parse_fasta(self):
        res = parse_fasta(
            [
                ">Seq1 abc def\n",
                "GGCTGCTATCAG\n",
                "CTAGCATCGTCGCATCGAC\n",
                ">Seq2\n",
                "ACGCTAGCTGCAAAA\n",
            ]
        )
        self.assertEqual(next(res), ("Seq1 abc def", "GGCTGCTATCAGCTAGCATCGTCGCATCGAC"))
        self.assertEqual(next(res), ("Seq2", "ACGCTAGCTGCAAAA"))
        self.assertRaises(StopIteration, next, res)

    def test_parse_empty_fasta(self):
        res = parse_fasta([])
        list_res = list(res)
        self.assertEqual(list_res, [])

    def test_load_fasta(self):
        f = tempfile.NamedTemporaryFile(mode="wt", encoding="utf-8")
        f.write(">Myseq asdf\n" "GGCTAAGGCCT\n" ">2ndseq *@#\n" "CCCGG\n")
        f.seek(0)
        self.assertEqual(
            load_fasta(f.name), {"Myseq": "GGCTAAGGCCT", "2ndseq": "CCCGG"}
        )

    def test_write_fasta(self):
        f = tempfile.NamedTemporaryFile(mode="w+t", encoding="utf-8")
        seqs = [("a", "CCGGT"), ("b", "TTTTTTTTT")]
        write_fasta(f, seqs)
        f.seek(0)
        self.assertEqual(f.read(), ">a\nCCGGT\n>b\nTTTTTTTTT\n")


class ResultsTests(unittest.TestCase):
    def test_parse_results(self):
        results = [
            "query_id\tspecies\ttypestrain_id\tregion_mismatches\tregion_positions\tprobability_incompatible\n",
            "Seq1\tA\tB\t1\t2\t0.5\n",
            "Seq1\tC\tD\t3\t1500\t-2.062794379753541e-12\n",
            "Seq2\tNA\tNA\tNA\tNA\tNA\n",
        ]

        with tempfile.NamedTemporaryFile(mode="w+t", encoding="utf-8") as f:
            f.writelines(results)
            f.seek(0)
            content = list(parse_results((l for l in f.readlines())))


if __name__ == "__main__":
    unittest.main()
