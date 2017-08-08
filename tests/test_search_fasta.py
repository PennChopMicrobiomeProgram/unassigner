from collections import namedtuple
import os.path
import tempfile
import unittest

from unassign.download import make_blast_db
from unassign.search_fasta import (
    FastaAligner, FastaAlignment,
    )

DATA_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "data")
GG_FP = os.path.join(DATA_DIR, "gg10.fasta")

class FastaAlignerTests(unittest.TestCase):
    def test_search_species(self):
        aligner = FastaAligner(GG_FP)
        seqs = [
            ("a", "CTTGCTCTCGGGTGACGAGCGGCGGACGGGTGAGTAAT"),
            ("b", "GCGTGGCGAACGGCTGACGAACACGTGG"),
            ]
        hits = aligner.search_species(seqs)
        observed = [(hit.query_id, hit.subject_id) for hit in hits]
        expected = [('b', '5'), ('b', '2'), ('b', '10'), ('b', '7'), ('b', '1')]
        self.assertEqual(observed, expected)


