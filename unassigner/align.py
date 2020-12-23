import os.path
import subprocess
import tempfile

from unassigner.parse import write_fasta, load_fasta, parse_fasta
from unassigner.alignment import AlignedPair

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

    def search(
            self, seqs, input_fp=None, output_fp=None, **kwargs):
        """Write seqs to input file, search, and parse output
        """
        if input_fp is None:
            infile = tempfile.NamedTemporaryFile(mode="w+t", encoding="utf-8")
            write_fasta(infile, seqs)
            infile.seek(0)
            input_fp = infile.name
        else:
            with open(input_fp, "w") as f:
                write_fasta(f, seqs)

        if output_fp is None:
            outfile = tempfile.NamedTemporaryFile()
            output_fp = outfile.name

        self._call(input_fp, output_fp, **kwargs)

        with open(output_fp) as f:
            for hit in self.parse(f):
                yield hit

    def _call(
            self, query_fp, output_fp, min_id = 0.5,
            maxaccepts = 1, threads=None, top_hits_only=False):
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
        return self._parse(f, self.convert_types, self.fields)

    @classmethod
    def _parse(self, f, convert_types=True, fields=DEFAULT_BLAST_FIELDS):
        """Parse a BLAST output file."""
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            vals = line.split("\t")
            res = dict(zip(fields, vals))
            if convert_types:
                for field in fields:
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


class HitExtender:
    def __init__(self, query_seqs, ref_seqs):
        self.query_seqs = dict(query_seqs)
        self.ref_seqs = dict(ref_seqs)

    def extend_hit(self, hit):
        # Handle the simple case where the local alignment covers both
        # sequences completely
        if self._is_global(hit):
            return AlignedPair(
                (hit['qseqid'], hit['qseq']),
                (hit['sseqid'], hit['sseq']))

        # We are going to need some repair or realignment.
        qseq = self.query_seqs[hit['qseqid']]
        assert(len(qseq) == hit['qlen'])
        sseq = self.ref_seqs[hit['sseqid']]
        assert(len(sseq) == hit['slen'])

        # We expect this always to be true for the vsearch aligner.
        # In the past, we had to re-do the alignment ourselves if it
        # didn't come back from the external program as a semi-global
        # alignment. I'm leaving this in place to double-check that we
        # have a semi-global alignment. If we go around adding endgaps
        # to an alignment that's not semi-global, we are sure to
        # generate some difficult-to-diagnose bugs.
        assert self._is_semiglobal(hit)

        qleft, sleft = self._add_endgaps_left(hit, qseq, sseq)
        qright, sright = self._add_endgaps_right(hit, qseq, sseq)
        aligned_qseq = qleft + hit['qseq'] + qright
        aligned_sseq = sleft + hit['sseq'] + sright
        return AlignedPair(
                (hit['qseqid'], aligned_qseq),
                (hit['sseqid'], aligned_sseq))

    @staticmethod
    def _is_global(hit):
        return (
            (hit['qstart'] == 1) and \
            (hit['sstart'] == 1) and \
            (hit['qend'] == hit['qlen']) and \
            (hit['send'] == hit['slen']))

    @staticmethod
    def _is_semiglobal(hit):
        more_to_the_left = (hit['qstart'] > 1) and \
                           (hit['sstart'] > 1)
        more_to_the_right = (hit['qend'] < hit['qlen']) and \
                            (hit['send'] < hit['slen'])
        return not (more_to_the_left or more_to_the_right)

    @staticmethod
    def _add_endgaps_left(hit, qseq, sseq):
        # No repair needed
        if (hit['qstart'] == 1) and (hit['sstart'] == 1):
            return ("", "")
        # Query hanging off to the left
        if (hit['qstart'] > 1) and (hit['sstart'] == 1):
            endgap_len = hit['qstart'] - 1
            return (qseq[:endgap_len], "-" * endgap_len)
        # Subject hanging off to the left
        if (hit['qstart'] == 1) and (hit['sstart'] > 1):
            endgap_len = hit['sstart'] - 1
            return ("-" * endgap_len, sseq[:endgap_len])
        # Anything not meeting these conditions is bad
        if (hit['qstart'] > 1) and (hit['sstart'] > 1):
            raise ValueError("Unaligned sequence on left")
        raise ValueError("Query or subject start position less than 1")

    @staticmethod
    def _add_endgaps_right(hit, qseq, sseq):
        # No repair needed
        if (hit['qend'] == hit['qlen']) and (hit['send'] == hit['slen']):
            return ("", "")
        # Query hanging off to the right
        if (hit['qend'] < hit['qlen']) and (hit['send'] == hit['slen']):
            endgap_len = hit['qlen'] - hit['qend']
            return (qseq[-endgap_len:], "-" * endgap_len)
        # Subject hanging off to the right
        if (hit['qend'] == hit['qlen']) and (hit['send'] < hit['slen']):
            endgap_len = hit['slen'] - hit['send']
            return ("-" * endgap_len, sseq[-endgap_len:])
        # Anything not meeting these conditions is bad
        if (hit['qend'] < hit['qlen']) and (hit['send'] < hit['qlen']):
            raise ValueError("Unaligned sequence on right")
        raise ValueError("Query or subject end position greater than length")

    def _get_subject_seq(self, subject_id):
        subject_outfile = tempfile.NamedTemporaryFile()
        subject_outfile_fp = subject_outfile.name
        args = ["blastdbcmd",
                "-db", self.db,
                "-entry", subject_id,
                "-out", subject_outfile_fp
        ]
        subprocess.check_call(args)
        with open(subject_outfile_fp) as f:
            return list(parse_fasta(f, trim_desc=True))[0][1]
