import io
import os.path
import unittest
import shutil
import tempfile

from unassigner.ani import (
    RefseqAssembly, Refseq16SDatabase,
)

DATA_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "data")
SUMMARY_FP = os.path.join(DATA_DIR, "assembly_summary.txt")


class RefseqAssemblyTests(unittest.TestCase):
    def setUp(self):
        self.a = RefseqAssembly(
            assembly_accession = "GCF_000010525.1",
            ftp_path = "ftp://ncbi.gov/abcd",
            asm_name = "ASM1052v1")

    def test_parse_summary(self):
        with open(SUMMARY_FP) as f:
            assemblies = list(RefseqAssembly.parse_summary(f))
            self.assertEqual(assemblies[0].accession, "GCF_000010525.1")
            self.assertEqual(assemblies[0].ftp_path, (
                "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/"
                "GCF/000/010/525/GCF_000010525.1_ASM1052v1"))

    def test_rna_url(self):
        self.assertEqual(
            self.a.rna_url,
            "https://ncbi.gov/abcd/abcd_rna_from_genomic.fna.gz")

    def test_genome_url(self):
        self.assertEqual(
            self.a.genome_url,
            "https://ncbi.gov/abcd/abcd_genomic.fna.gz")


class MockAssembly:
    def __init__(self, acc, ssu_seqs):
        self.accession = acc
        self.ssu_seqs = ssu_seqs


class Refseq16SDatabaseTests(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.db = Refseq16SDatabase()
        self.db.ssu_accession_fp = os.path.join(self.dir, "accessions.txt")
        self.db.ssu_fasta_fp = os.path.join(self.dir, "seqs.fasta")

        a_seqs = [("s1", "TCCG"), ("s1_duplicate", "TCCG"), ("s2", "TTTT")]
        self.a = MockAssembly("assembly1.1", a_seqs)

        self.bs = MockAssembly("bs", [("bs16S", bs_seq)])
        self.bv = MockAssembly("bv", [("bv16S", bv_seq)])
        self.sp = MockAssembly("sp", [("sp16S", sp_seq)])

    def tearDown(self):
        shutil.rmtree(self.dir)

    def test_add_assembly(self):
        db = Refseq16SDatabase()
        db.add_assembly(self.a)
        self.assertEqual(db.seqs, {"s1": "TCCG", "s2": "TTTT"})
        self.assertEqual(db.assemblies, {"s1": self.a, "s2": self.a})
        self.assertEqual(db.seqids_by_assembly["assembly1.1"], ["s1", "s2"])

    def test_load_from_objects(self):
        assemblies = {"bv": self.bv, "sps": self.sp}
        was_loaded_from_file = self.db.load(assemblies)

        self.assertFalse(was_loaded_from_file)
        self.assertEqual(self.db.seqs["bv16S"], bv_seq)
        self.assertEqual(self.db.seqs["sp16S"], sp_seq)
        self.assertEqual(self.db.assemblies, {
            "bv16S": self.bv,
            "sp16S": self.sp,
        })
        self.assertEqual(self.db.seqids_by_assembly["bv"], ["bv16S"])
        self.assertEqual(self.db.seqids_by_assembly["sp"], ["sp16S"])

        with open(self.db.ssu_accession_fp) as f:
            observed_accession_txt = f.read()
        self.assertEqual(observed_accession_txt, ACCESSION_TXT)
        with open(self.db.ssu_fasta_fp) as f:
            observed_seqs_txt = f.read()
        self.assertEqual(observed_seqs_txt, SEQS_TXT)

    def test_load_from_file(self):
        with open(self.db.ssu_accession_fp, "wt") as f:
            f.write(ACCESSION_TXT)
        with open(self.db.ssu_fasta_fp, "wt") as f:
            f.write(SEQS_TXT)

        # Ensure we are loading from files, not objects
        self.bv.ssu_seqs = None
        self.sp.ssu_seqs = None

        assemblies = {"bv": self.bv, "sp": self.sp}
        was_loaded_from_file = self.db.load(assemblies)

        self.assertTrue(was_loaded_from_file)
        self.assertEqual(self.db.seqs["bv16S"], bv_seq)
        self.assertEqual(self.db.seqs["sp16S"], sp_seq)
        self.assertEqual(self.db.assemblies, {
            "bv16S": self.bv,
            "sp16S": self.sp,
        })

    def test_search_one(self):
        self.db.load({
            "bv": self.bv, "sp": self.sp, "bs": self.bs,
        })
        all_results = list(self.db.search_at_pctid("bv16S", 96.5))
        self.assertEqual(len(all_results), 1)

        r = all_results[0]
        self.assertEqual(r.query, self.bv)
        self.assertEqual(r.subject, self.bs)
        self.assertEqual(r.pctid, 96.5)
        self.assertEqual(r.query_seqid, "bv16S")
        self.assertEqual(r.subject_seqid, "bs16S")

    def test_save_accessions(self):
        db = Refseq16SDatabase()
        db.add_assembly(self.a)
        f = io.StringIO()
        db._save_accessions(f)
        self.assertEqual(f.getvalue(), "s1\tassembly1.1\ns2\tassembly1.1\n")

    def test_load_accessions(self):
        accession_lines = [
            "s1\tassembly1.1\n",
            "s2\tassembly1.1\n",
        ]
        db = Refseq16SDatabase()
        db._load_accessions(accession_lines, {self.a.accession: self.a})
        self.assertEqual(db.assemblies, {"s1": self.a, "s2": self.a})
        self.assertEqual(db.seqids_by_assembly["assembly1.1"], ["s1", "s2"])
        self.assertEqual(db.seqs, {})

    def test_save_seqs(self):
        db = Refseq16SDatabase()
        db.add_assembly(self.a)
        f = io.StringIO()
        db._save_seqs(f)
        self.assertEqual(f.getvalue(), ">s1\nTCCG\n>s2\nTTTT\n")

    def test_load_seqs(self):
        seq_lines = [">s1\n", "TCCG\n", ">s2\n", "TTTT\n"]
        db = Refseq16SDatabase()
        db._load_seqs(seq_lines)
        self.assertEqual(db.assemblies, {})
        self.assertEqual(db.seqids_by_assembly["assembly1.1"], [])
        self.assertEqual(db.seqs, {"s1": "TCCG", "s2": "TTTT"})

bs_seq = "AGAGTTTGATCCTGGCTCAGGATGAACGCTAGCTACAGGCTTAACACATGCAAGTCGAGGGGCAGCATGGTCTTAGCTTGCTAAGGCTGATGGCGACCGGCGCACGGGTGAGTAACACGTATCCAACCTGCCGTCTACTCTTGGCCAGCCTTCTGAAAGGAAGATTAATCCAGGATGGGATCATGAGTTCACATGTCCGCATGATTAAAGGTATTTTCCGGTAGACGATGGGGATGCGTTCCATTAGATAGTAGGCGGGGTAACGGCCCACCTAGTCAACGATGGATAGGGGTTCTGAGAGGAAGGTCCCCCACATTGGAACTGAGACACGGTCCAAACTCCTACGGGAGGCAGCAGTGAGGAATATTGGTCAATGGGCGATGGCCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATAAAGGAATAAAGTCGGGTATGCATACCCGTTTGCATGTACTTTATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGATGGATGTTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGATACTGGATGTCTTGAGTGCAGTTGAGGCAGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCCTGCTAAGCTGCAACTGACATTGAGGCTCGAAAGTGTGGGTATCAAACAGGATTAGATACCCTGGTAGTCCACACGGTAAACGATGAATACTCGCTGTTTGCGATATACGGCAAGCGGCCAAGCGAAAGCGTTAAGTATTCCACCTGGGGAGTACGCCGGCAACGGTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGAGGAACATGTGGTTTAATTCGATGATACGCGAGGAACCTTACCCGGGCTTAAATTGCACTCGAATGATCCGGAAACGGTTCAGCTAGCAATAGCGAGTGTGAAGGTGCTGCATGGTTGTCGTCAGCTCGTGCCGTGAGGTGTCGGCTTAAGTGCCATAACGAGCGCAACCCTTGTTGTCAGTTACTAACAGGTGATGCTGAGGACTCTGACAAGACTGCCATCGTAAGATGTGAGGAAGGTGGGGATGACGTCAAATCAGCACGGCCCTTACGTCCGGGGCTACACACGTGTTACAATGGGGGGTACAGAGGGCCGCTACCACGCGAGTGGATGCCAATCCCTAAAACCCCTCTCAGTTCGGACTGGAGTCTGCAACCCGACTCCACGAAGCTGGATTCGCTAGTAATCGCGCATCAGCCACGGCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCAAGCCATGGGAGCCGGGGGTACCTGAAGTGCGTAACCGCGAGGATCGCCCTAGGGTAAAACTGGTGACTGGGGCTAAGTCTAACCAAGGTAACC"

bv_seq = "GTTTGATCCTGGCTCAGGATGAACGCTAGCTACAGGCTTAACACATGCAAGTCGAGGGGCAGCATGGTCTTAGCTTGCTAAGGCCGATGGCGACCGGCGCACGGGTGAGTAACACGTATCCAACCTGCCGTCTACTCTTGGACAGCCTTCTGAAAGGAAGATTAATACAAGATGGCATCATGAGTCCGCATGTTCACATGATTAAAGGTATTCCGGTAGACGATGGGGATGCGTTCCATTAGATAGTAGGCGGGGTAACGGCCCACCTAGTCTTCGATGGATAGGGTTCTGAGAGGAAGTCCCCCACATTGGAACTGAGACACGGTCCAAACTCCTACGGGAGGCAGCAGTGAGGAATATTGGTCAATGGGCGAGAGCCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATAAAGGAATAAAGTCGGGTATGGATACCCGTTTGCATGTACTTTATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGATGGATGTTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGATACTGGATATCTTGAGTGCAGTTGAGGCAGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCCTGCTAAGCTGCAACTGACATTGAGGCTCGAAAGTGTGGGTATCAAACAGGATTAGATACCCTGGTAGTCCACACGGTAAACGATGAATACTCGCTGTTTGCGATATACGGCAAGCGGCCAAGCGAAAGCGTTAAGTATTCCACCTGGGGAGTACGCCGGCAACGGTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGAGGAACATGTGGTTTAATTCGATGATACGCGAGGAACCTTACCCGGGCTTAAATTGCAGATGAATTACGGTGAAAGCCGTAAGCCGCAAGGCATCTGTGAAGGTGCTGCATGGTTGTCGTCAGCTCGTGCCGTGAGGTGTCGGCTTAAGTGCCATAACGAGCGCAACCCTTGTTGTCAGTTACTAACAGGTTTTGCTGAGGACTCTGACAAGACTGCCATCGTAAGATGTGAGGAAGGTGGGGATGACGTCAAATCAGCACGGCCCTTACGTCCGGGGCTACACACGTGTTACAATGGGGGGTACAGAGGGCCGCTACCACGCGAGTGATGCCAATCCCCAAAACCTCTCTCAGTTCGGACTGGAGTCTGCAACCCGACTCCACGAAGCTGGATTCGCTAGTAATCGCGCATCAGCCACGGCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCAAGCCATGGGAGCCGGGGGGTACCTGAAGTGCGTAACCGCGAGGAGCGCCCTAGGGTAAAACTGGTGACTGGGGCTAAGTCGTAACA"

sp_seq = "GAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCTGTAAACGATGAGTGCTAGGTGTTAGACCCTTTCCGGGGTTTAGTGCCGTAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGACCGCAAGGTTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCCTCTGACCGCTCTAGAGATAGAGTTTTCCTTCGGGACAGAGGTGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCCTATTGTTAGTTGCCATCATTCAGTTGGGCACTCTAGCGAGACTGCCGGTAATAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGACCTGGGCTACACACGTGCTACAATGGCTGGTACAACGAGTCGCAAGCCGGTGACGGCAAGCTAATCTCTTAAAGCCAGTCTCAGTTCGGATTGTAGGCTGCAACTCGCCTACATGAAGTCGGAATCGCTAGTAATCGCGGATCAGCACGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAGAGTTTGTAACACCCGAAGTCGGTGAGGTAACCGTAAGGAGCCAGCCGCCTAAGGTGGGATAGATGATTGGGGTGAAGTCGTAACAAGGTCAGCCGTTTGGGAGA"

ACCESSION_TXT = """\
bv16S\tbv
sp16S\tsp
"""

SEQS_TXT = ">bv16S\n{0}\n>sp16S\n{1}\n".format(bv_seq, sp_seq)
