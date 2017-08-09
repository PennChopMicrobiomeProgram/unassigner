import subprocess

from unassign.alignment import Aligner


class FastaAlignment(object):
    def __init__(self, hit):
        self._hit = hit
        self.query_id = hit["query_id"]
        self.subject_id = hit["subject_id"]


class FastaAligner(Aligner):
    alignment_cls = FastaAlignment
    executable = "glsearch36"
    field_types = [
        str, str, float, int,
        int, int, int, int, int, int,
        float, float, str,
    ]
    field_names = [
        "query_id", "subject_id", "pct_identity", "alignment_length",
        "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end",
        "evalue", "bit_score", "BTOP",
    ]

    def _call(self, query_fp, library_fp, output_fp):
        args = [
            "glsearch36", "-n",
            "-T", str(self.num_cpus),
            "-b", str(self.species_max_hits),
            "-m", "8B",
            query_fp, library_fp,
            ]
        # Noticed that query seqs were missing with -O flag
        # Open file and redirect stdout instead
        with open(output_fp, "w") as f:
            subprocess.check_call(args, stdout=f)

