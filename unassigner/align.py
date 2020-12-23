import os.path
import subprocess
import tempfile

from unassigner.parse import write_fasta

DEFAULT_BLAST_FIELDS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "qlen", "slen", "qseq", "sseq",
]

BLAST_FIELD_TYPES = {
    "qseqid": str,
    "sseqid": str,
    "pident": float,
    "length": int,
    "mismatch": int,
    "gapopen": int,
    "qstart": int,
    "qend": int,
    "sstart": int,
    "send": int,
    "qlen": int,
    "slen": int,
    "qseq": str,
    "sseq": str,
}

VSEARCH_TO_BLAST = {
    "query": "qseqid",
    "target": "sseqid",
    "id2": "pident",
    "alnlen": "length",
    "mism": "mismatch",
    "gaps": "gapopen",
    "qilo": "qstart",
    "qihi": "qend",
    "tilo": "sstart",
    "tihi": "send",
    "qs": "qlen",
    "ts": "slen",
    "qrow": "qseq",
    "trow": "sseq",
}

BLAST_TO_VSEARCH = {b: v for v, b in VSEARCH_TO_BLAST.items()}


class VsearchAligner:
    def __init__(self, ref_seqs_fp):
        self.ref_seqs_fp = ref_seqs_fp
        self.fields = DEFAULT_BLAST_FIELDS
        self.convert_types = True
        self.stderr = subprocess.DEVNULL

    def search(self, seqs, input_fp=None, output_fp=None, **kwargs):
        """Search seqs and return hits"""
        if input_fp is None:
            infile = tempfile.NamedTemporaryFile(
                suffix=".fasta", mode="w+t", encoding="utf-8")
            write_fasta(infile, seqs)
            infile.seek(0)
            input_fp = infile.name
        else:
            with open(input_fp, "w") as f:
                write_fasta(f, seqs)

        if output_fp is None:
            outfile = tempfile.NamedTemporaryFile(
                suffix=".txt", mode="wt")
            output_fp = outfile.name

        self._call(input_fp, output_fp, **kwargs)

        with open(output_fp) as f:
            for hit in self.parse(f):
                yield hit

    def _call(
            self, query_fp, output_fp, min_id = 0.5,
            maxaccepts = 1, threads=None, top_hits_only=False):
        """Call vsearch to generate output file"""
        self.make_reference_udb()
        id_arg = "{:.3f}".format(min_id)
        userfields_arg = "+".join(BLAST_TO_VSEARCH[f] for f in self.fields)
        maxaccepts_arg = "{:d}".format(maxaccepts)
        args = [
            "vsearch", "--usearch_global", query_fp,
            "--db", self.ref_seqs_udb_fp,
            "--iddef", "2",
            "--id", id_arg,
            "--userout", output_fp,
            "--userfields", userfields_arg,
            "--maxaccepts", maxaccepts_arg,
        ]
        if threads is not None:
            threads_arg = "{:d}".format(threads)
            args.extend(["--threads", threads_arg])
        if top_hits_only:
            args.append("--top_hits_only")
        subprocess.check_call(args, stderr=self.stderr)

    def parse(self, f):
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            vals = line.split("\t")
            res = dict(zip(self.fields, vals))
            if self.convert_types:
                for field in self.fields:
                    fcn = BLAST_FIELD_TYPES[field]
                    res[field] = fcn(res[field])
            yield res

    @property
    def ref_seqs_udb_fp(self):
        base_fp, _ = os.path.splitext(self.ref_seqs_fp)
        return base_fp + ".udb"

    def make_reference_udb(self):
        if os.path.exists(self.ref_seqs_udb_fp):
            return
        args = [
            "vsearch",
            "--makeudb_usearch", self.ref_seqs_fp,
            "--output", self.ref_seqs_udb_fp,
        ]
        return subprocess.check_call(args, stderr=self.stderr)
