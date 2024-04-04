import collections
import os
import io
import tempfile
import unittest

from unassigner.mismatch_db import (
    MismatchLocationApp,
    main,
    group_by_n,
    MutableMismatchDb,
)

DATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")


class MismatchLocationAppTests(unittest.TestCase):
    def setUp(self):
        self.oral_species_fp = os.path.join(DATA_DIR, "oral_species.fasta")
        self.oral_reference_fp = os.path.join(DATA_DIR, "oral_refs.fasta")
        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        udb_fp = os.path.join(DATA_DIR, "oral_refs.udb")
        if os.path.exists(udb_fp):
            os.remove(udb_fp)

    def test_mismatch_location_app(self):
        output_file = io.StringIO()
        with open(self.oral_species_fp) as species_file:
            app = MismatchLocationApp(species_file, self.oral_reference_fp, output_file)
            app.run()
        self.assertEqual(output_file.getvalue(), oral_mismatches)

    def test_mismatch_command_line(self):
        output_fp = os.path.join(self.dir, "mismatches.txt")
        main(
            [
                self.oral_species_fp,
                self.oral_reference_fp,
                output_fp,
            ]
        )
        with open(output_fp) as f:
            output_txt = f.read()
        self.assertEqual(output_txt, oral_mismatches)

    def test_mismatch_location_app_batch_size(self):
        output_file = io.StringIO()
        with open(self.oral_species_fp) as species_file:
            app = MismatchLocationApp(
                species_file, self.oral_reference_fp, output_file, batch_size=2
            )
            app.run()
        self.assertEqual(output_file.getvalue(), oral_mismatches)


class MismatchDbFunctionTests(unittest.TestCase):
    def test_group_by_n(self):
        self.assertEqual(
            list(group_by_n([1, 2, 3, 4, 5, 6, 7, 8], 3)),
            [[1, 2, 3], [4, 5, 6], [7, 8]],
        )


oral_mismatches = """\
AF003928	4446824	628	956	1033
AF003928	4439903	135	357	1076
AF003928	4367196	763	1130	1386
AF003928	1057832	305	369	821
AF003928	983604	655	703	1130
AF003928	894424	147	206	601
AF003928	735302	147	206	601
AF003928	691335	305	369	821
AF003928	1103279	206	211	689
AF003928	1101991	206	221	1130
"""

oral_mismatch_data = {
    "AF003928": [
        ("4446824", [628, 956, 1033]),
        ("4439903", [135, 357, 1076]),
        ("4367196", [763, 1130, 1386]),
        ("1057832", [305, 369, 821]),
        ("983604", [655, 703, 1130]),
        ("894424", [147, 206, 601]),
        ("735302", [147, 206, 601]),
        ("691335", [305, 369, 821]),
        ("1103279", [206, 211, 689]),
        ("1101991", [206, 221, 1130]),
    ]
}
