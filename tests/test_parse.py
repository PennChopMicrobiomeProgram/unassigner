import unittest

from unassign.parse import parse_fasta

class FastaTests(unittest.TestCase):
    def test_parse_fasta(self):
        res = parse_fasta([
            ">Seq1 abc def\n",
            "GGCTGCTATCAG\n",
            "CTAGCATCGTCGCATCGAC\n",
            ">Seq2\n",
            "ACGCTAGCTGCAAAA\n",
            ])
        self.assertEqual(next(res), (
            "Seq1 abc def", "GGCTGCTATCAGCTAGCATCGTCGCATCGAC"))
        self.assertEqual(next(res), ("Seq2", "ACGCTAGCTGCAAAA"))
        self.assertRaises(StopIteration, next, res)

if __name__ == "__main__":
    unittest.main()
