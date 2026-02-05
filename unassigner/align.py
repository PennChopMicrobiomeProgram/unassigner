import abc
import os.path
import subprocess
import tempfile
from Bio.Align import PairwiseAligner

from unassigner.parse import write_fasta, parse_fasta
from unassigner.alignment import AlignedPair

BLAST_FMT = (
    "qseqid sseqid pident length mismatch gapopen "
    "qstart qend sstart send qlen slen qseq sseq"
)
BLAST_FIELDS = BLAST_FMT.split()
BLAST_FIELD_TYPES = [
    str,
    str,
    float,
    int,
    int,
    int,
    int,
    int,
    int,
    int,
    int,
    int,
    str,
    str,
]


def _hit_stats(aligned_qseq, aligned_sseq):
    matches = 0
    mismatches = 0
    gapopen = 0
    in_gap = False
    for qbase, sbase in zip(aligned_qseq, aligned_sseq):
        if qbase == "-" or sbase == "-":
            if not in_gap:
                gapopen += 1
                in_gap = True
            continue
        in_gap = False
        if qbase == sbase:
            matches += 1
        else:
            mismatches += 1
    alnlen = len(aligned_qseq)
    pident = (matches / alnlen) * 100 if alnlen else 0.0

    q_positions = [i for i, base in enumerate(aligned_qseq) if base != "-"]
    s_positions = [i for i, base in enumerate(aligned_sseq) if base != "-"]
    qstart = q_positions[0] + 1
    qend = q_positions[-1] + 1
    sstart = s_positions[0] + 1
    send = s_positions[-1] + 1
    qlen = len(q_positions)
    slen = len(s_positions)

    return {
        "pident": pident,
        "length": alnlen,
        "mismatch": mismatches,
        "gapopen": gapopen,
        "qstart": qstart,
        "qend": qend,
        "sstart": sstart,
        "send": send,
        "qlen": qlen,
        "slen": slen,
    }


class Aligner(abc.ABC):
    def __init__(self, ref_seqs_fp):
        self.ref_seqs_fp = ref_seqs_fp

    def search(self, seqs, input_fp=None, output_fp=None, **kwargs):
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

        self._call(input_fp, self.ref_seqs_fp, output_fp, **kwargs)

        with open(output_fp) as f:
            for hit in self._parse(f):
                yield hit

    @classmethod
    def _parse(self, f, convert_types=True):
        """Parse a BLAST output file."""
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            vals = line.split("\t")
            if convert_types:
                vals = [fn(v) for fn, v in zip(BLAST_FIELD_TYPES, vals)]
            yield dict(zip(BLAST_FIELDS, vals))


class VsearchAligner(Aligner):
    @property
    def ref_seqs_udb_fp(self):
        base_fp, _ = os.path.splitext(self.ref_seqs_fp)
        return base_fp + ".udb"

    def make_reference_udb(self):
        if os.path.exists(self.ref_seqs_udb_fp):
            return True
        args = [
            "vsearch",
            "--makeudb_usearch",
            self.ref_seqs_fp,
            "--output",
            self.ref_seqs_udb_fp,
        ]
        try:
            subprocess.check_call(args)
            return True
        except FileNotFoundError:
            return False

    @staticmethod
    def _index(fasta_fp):
        pass

    def _call(self, query_fp, database_fp, output_fp, min_id=0.5, **kwargs):
        """Call the VSEARCH program.

        --id is a required argument to run the program. We expose the
        argument here as min_id.
        """
        if not os.path.exists(self.ref_seqs_udb_fp):
            has_vsearch = self.make_reference_udb()
            if not has_vsearch:
                self._fallback_call(query_fp, output_fp, min_id=min_id, **kwargs)
                return
        args = [
            "vsearch",
            "--usearch_global",
            query_fp,
            "--db",
            self.ref_seqs_udb_fp,
            "--iddef",
            "2",
            "--id",
            str(min_id),
            "--userout",
            output_fp,
            "--userfields",
            "query+target+id2+alnlen+mism+gaps+qilo+qihi+tilo+tihi+qs+ts+qrow+trow",
        ]
        for arg, val in kwargs.items():
            arg = "--" + arg
            if val is None:
                args.append(arg)
            else:
                args += [arg, str(val)]
        try:
            subprocess.check_call(args, stderr=subprocess.DEVNULL)
        except FileNotFoundError:
            self._fallback_call(query_fp, output_fp, min_id=min_id, **kwargs)

    def _fallback_call(self, query_fp, output_fp, min_id=0.5, **kwargs):
        top_hits_only = "top_hits_only" in kwargs or kwargs.get("maxaccepts") == 1
        with open(query_fp) as qf:
            query_recs = list(parse_fasta(qf, trim_desc=True))
        with open(self.ref_seqs_fp) as sf:
            subject_recs = list(parse_fasta(sf, trim_desc=True))

        with open(output_fp, "w") as out:
            for query_id, query_seq in query_recs:
                query_hits = []
                for subject_id, subject_seq in subject_recs:
                    aligned_qseq, aligned_sseq = align_semiglobal(
                        query_seq, subject_seq
                    )
                    stats = _hit_stats(aligned_qseq, aligned_sseq)
                    if (stats["pident"] / 100) < min_id:
                        continue
                    hit = [
                        query_id,
                        subject_id,
                        stats["pident"],
                        stats["length"],
                        stats["mismatch"],
                        stats["gapopen"],
                        stats["qstart"],
                        stats["qend"],
                        stats["sstart"],
                        stats["send"],
                        stats["qlen"],
                        stats["slen"],
                        aligned_qseq,
                        aligned_sseq,
                    ]
                    query_hits.append(hit)

                query_hits.sort(key=lambda h: (h[2], h[3]), reverse=True)
                if top_hits_only and query_hits:
                    query_hits = [query_hits[0]]
                for hit in query_hits:
                    out.write("\t".join(str(v) for v in hit))
                    out.write("\n")


class BlastAligner(Aligner):
    @staticmethod
    def _index(fasta_fp):
        return subprocess.check_call(
            [
                "makeblastdb",
                "-dbtype",
                "nucl",
                "-in",
                fasta_fp,
            ],
            stdout=subprocess.DEVNULL,
        )

    def _call(self, query_fp, database_fp, output_fp, **kwargs):
        """Call the BLAST program."""
        args = [
            "blastn",
            "-evalue",
            "1e-5",
            "-outfmt",
            "6 " + BLAST_FMT,
        ]
        for arg, val in kwargs.items():
            arg = "-" + arg
            if val is None:
                args.append(arg)
            else:
                args += [arg, str(val)]
        args += [
            "-query",
            query_fp,
            "-db",
            database_fp,
            "-out",
            output_fp,
        ]
        subprocess.check_call(args)


class HitExtender:
    def __init__(self, query_seqs, ref_seqs):
        self.query_seqs = dict(query_seqs)
        self.ref_seqs = dict(ref_seqs)

    def extend_hit(self, hit):
        # Handle the simple case where the local alignment covers both
        # sequences completely
        if self._is_global(hit):
            return AlignedPair(
                (hit["qseqid"], hit["qseq"]), (hit["sseqid"], hit["sseq"])
            )

        # We are going to need some repair or realignment.
        qseq = self.query_seqs[hit["qseqid"]]
        assert len(qseq) == hit["qlen"]
        sseq = self.ref_seqs[hit["sseqid"]]
        assert len(sseq) == hit["slen"]

        if self._needs_realignment(hit):
            aligned_qseq, aligned_sseq = align_semiglobal(qseq, sseq)
            return AlignedPair(
                (hit["qseqid"], aligned_qseq), (hit["sseqid"], aligned_sseq)
            )

        qleft, sleft = self._add_endgaps_left(hit, qseq, sseq)
        qright, sright = self._add_endgaps_right(hit, qseq, sseq)
        aligned_qseq = qleft + hit["qseq"] + qright
        aligned_sseq = sleft + hit["sseq"] + sright
        return AlignedPair((hit["qseqid"], aligned_qseq), (hit["sseqid"], aligned_sseq))

    @staticmethod
    def _is_global(hit):
        return (
            (hit["qstart"] == 1)
            and (hit["sstart"] == 1)
            and (hit["qend"] == hit["qlen"])
            and (hit["send"] == hit["slen"])
        )

    @staticmethod
    def _needs_realignment(hit):
        more_to_the_left = (hit["qstart"] > 1) and (hit["sstart"] > 1)
        more_to_the_right = (hit["qend"] < hit["qlen"]) and (hit["send"] < hit["slen"])
        return more_to_the_left or more_to_the_right

    @staticmethod
    def _add_endgaps_left(hit, qseq, sseq):
        # No repair needed
        if (hit["qstart"] == 1) and (hit["sstart"] == 1):
            return ("", "")
        # Query hanging off to the left
        if (hit["qstart"] > 1) and (hit["sstart"] == 1):
            endgap_len = hit["qstart"] - 1
            return (qseq[:endgap_len], "-" * endgap_len)
        # Subject hanging off to the left
        if (hit["qstart"] == 1) and (hit["sstart"] > 1):
            endgap_len = hit["sstart"] - 1
            return ("-" * endgap_len, sseq[:endgap_len])
        # Anything not meeting these conditions is bad
        if (hit["qstart"] > 1) and (hit["sstart"] > 1):
            raise ValueError("Unaligned sequence on left")
        raise ValueError("Query or subject start position less than 1")

    @staticmethod
    def _add_endgaps_right(hit, qseq, sseq):
        # No repair needed
        if (hit["qend"] == hit["qlen"]) and (hit["send"] == hit["slen"]):
            return ("", "")
        # Query hanging off to the right
        if (hit["qend"] < hit["qlen"]) and (hit["send"] == hit["slen"]):
            endgap_len = hit["qlen"] - hit["qend"]
            return (qseq[-endgap_len:], "-" * endgap_len)
        # Subject hanging off to the right
        if (hit["qend"] == hit["qlen"]) and (hit["send"] < hit["slen"]):
            endgap_len = hit["slen"] - hit["send"]
            return ("-" * endgap_len, sseq[-endgap_len:])
        # Anything not meeting these conditions is bad
        if (hit["qend"] < hit["qlen"]) and (hit["send"] < hit["qlen"]):
            raise ValueError("Unaligned sequence on right")
        raise ValueError("Query or subject end position greater than length")

    def _get_subject_seq(self, subject_id):
        subject_outfile = tempfile.NamedTemporaryFile()
        subject_outfile_fp = subject_outfile.name
        args = [
            "blastdbcmd",
            "-db",
            self.db,
            "-entry",
            subject_id,
            "-out",
            subject_outfile_fp,
        ]
        subprocess.check_call(args)
        with open(subject_outfile_fp) as f:
            return list(parse_fasta(f, trim_desc=True))[0][1]


def align_semiglobal(qseq, sseq):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 5
    aligner.mismatch_score = -4
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.end_open_gap_score = 0
    aligner.end_extend_gap_score = 0
    alignments = aligner.align(sseq, qseq)
    alignment = alignments[0]

    return alignment[1], alignment[0]
