import subprocess
import tempfile

import unassign.parse

class FastaAlignment(object):
    def __init__(self, hit):
        self._hit = hit
        self.query_id = hit["query_id"]
        self.subject_id = hit["subject_id"]

class FastaAligner(object):
    alignment_cls = FastaAlignment
    executable = "glsearch36"

    def __init__(self, species_fp):
        self.species_fp = species_fp
        self.species_max_hits = 5
        self.species_input_fp = None
        self.species_output_fp = None
        self.num_cpus = 1

    def search_species(self, seqs):
        """Search species typestrains for match to query sequences."""
        return self._search(
            seqs, self.species_fp, self.species_max_hits,
            self.species_input_fp, self.species_output_fp)

    def _search(self, seqs, db, max_hits, input_fp, output_fp):
        if input_fp is None:
            infile = tempfile.NamedTemporaryFile(mode="w+t", encoding="utf-8")
            unassign.parse.write_fasta(infile, seqs)
            infile.seek(0)
            input_fp = infile.name
        else:
            with open(input_fp, "w") as f:
                unassign.parse.write_fasta(f, seqs)

        if output_fp is None:
            outfile = tempfile.NamedTemporaryFile()
            output_fp = outfile.name

        self._call(input_fp, db, output_fp, max_hits, self.num_cpus)
        return self._load(output_fp)

    @classmethod
    def _load(self, output_fp):
        """Load hits from an output file."""
        with open(output_fp) as f:
            hits = [self.alignment_cls(x) for x in self._parse(f)]
        return hits

    @classmethod
    def _parse(self, f):
        FASTA_FIELD_TYPES = [
            str, str, float, int,
            int, int, int, int, int, int,
            float, float, str,
        ]
        FASTA_FIELDS = [
            "query_id", "subject_id", "pct_identity", "alignment_length",
            "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end",
            "evalue", "bit_score", "BTOP",
        ]
        for line in f:
            line = line.strip()
            if line.startswith("#") or (not line):
                continue
            vals = line.split("\t")
            vals = [fn(v) for fn, v in zip(FASTA_FIELD_TYPES, vals)]
            yield dict(zip(FASTA_FIELDS, vals))

    def _call(self, query_fp, library_fp, output_fp, max_target_seqs, num_cpus):
        args = [
            "glsearch36", "-n",
            "-T", str(num_cpus),
            "-b", str(max_target_seqs),
            "-m", "8B",
            query_fp, library_fp,
            ]
        # Noticed that query seqs were missing with -O flag
        # Open file and redirect stdout instead
        with open(output_fp, "w") as f:
            subprocess.check_call(args, stdout=f)

