import unittest

from unassigner.download import parse_desc, parse_fasta


class DownloadTests(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_parsers(self):
        inputs = (
            ">lcl\|LTP_06_2022\|TraGuam2 [organism=Organism name] [strain=NBRC 103172 s[T]] [accession=AC123456]\n"
            "AUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAGCGGCAGCGGGGGAA\n"
            ">lcl\|LTP_06_2022\|TraOdont [organism=Organism name] [strain=Eant 3-9 l[T] r[T] s[T]] [accession=AC123456]\n"
            "AUGCAAGUCGAGCGGCAGCGGGGGAAAGCUUGCUUUCCCGCCGGCGAGCGGCGG\n"
            ">lcl\|LTP_06_2022\|EnbTribo [organism=Organism name] [strain=IG-V01 l[T]] [accession=AC123456]\n"
            "UCCAGAGUUUGAUCAUGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGC\n"
        )

        for desc, seq in parse_fasta(inputs):
            print(desc)
            self.assertEqual(desc[0], ">")
            self.assertEqual(set(list(seq)), set("A", "C", "G", "U"))

            accession, species_name = parse_desc(desc)
            self.assertEqual(accession, "AC123456")
            self.assertEqual(species_name, "Organism name")
