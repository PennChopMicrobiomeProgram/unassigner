import tempfile
import unittest

from unassigner.parse import parse_fasta, parse_desc, load_fasta, write_fasta


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

    def test_parse_desc(self):
        seqs = {
            ">lcl\\|LTP_06_2022\\|TraGuam2 [organism=Organism name] [strain=NBRC 103172 s[T]] [accession=AC123456]\n": "AUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAGCGGCAGCGGGGGAA\n",
            ">lcl\\|LTP_06_2022\\|TraOdont [organism=Organism name] [strain=Eant 3-9 l[T] r[T] s[T]] [accession=AC123456]\n": "AUGCAAGUCGAGCGGCAGCGGGGGAAAGCUUGCUUUCCCGCCGGCGAGCGGCGG\n",
            ">lcl\\|LTP_06_2022\\|EnbTribo [organism=Organism name] [strain=IG-V01 l[T]] [accession=AC123456]\n": "UCCAGAGUUUGAUCAUGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGC\n",
        }

        for desc, seq in seqs.items():
            self.assertEqual(desc[0], ">")
            self.assertEqual(set(list(seq.strip("\n"))), set(["A", "C", "G", "U"]))

            accession, species_name = parse_desc(desc)
            self.assertEqual(accession, "AC123456")
            self.assertEqual(species_name, "Organism name")


if __name__ == "__main__":
    unittest.main()
