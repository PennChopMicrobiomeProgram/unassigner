import abc
import subprocess
import tempfile
from Bio import pairwise2

from unassigner.parse import write_fasta, load_fasta, parse_fasta
from unassigner.alignment import AlignedPair

BLAST_FMT = (
    "qseqid sseqid pident length mismatch gapopen "
    "qstart qend sstart send qlen slen qseq sseq")
BLAST_FIELDS = BLAST_FMT.split()
BLAST_FIELD_TYPES = [
    str, str, float, int, int, int,
    int, int, int, int, int, int, str, str]

# TODO: Need a BlastFilter to collect multiple hits and get the best one

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
    def _parse(self, f):
        """Parse a BLAST output file."""
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            vals = line.split("\t")
            vals = [fn(v) for fn, v in zip(BLAST_FIELD_TYPES, vals)]
            yield dict(zip(BLAST_FIELDS, vals))


class VsearchAligner(Aligner):
    @staticmethod
    def _index(fasta_fp):
        pass

    def _call(self, query_fp, database_fp, output_fp, min_id = 0.5, **kwargs):
        """Call the VSEARCH program.

        --id is a required argument to run the program. We expose the
        argument here as min_id.
        """
        args = [
            "vsearch", "--usearch_global", query_fp,
            "--db", database_fp,
            "--iddef", "2",
            "--id", str(min_id),
            "--userout", output_fp,
            "--userfields",
            "query+target+id2+alnlen+mism+gaps+qilo+qihi+tilo+tihi+qs+ts+qrow+trow"
            ]
        for arg, val in kwargs.items():
            arg = "--" + arg
            if val is None:
                args.append(arg)
            else:
                args += [arg, str(val)]
        subprocess.check_call(args, stderr=subprocess.DEVNULL)


class BlastAligner(Aligner):
    @staticmethod
    def _index(fasta_fp):
        return subprocess.check_call([
            "makeblastdb",
            "-dbtype", "nucl",
            "-in", fasta_fp,
            ], stdout=subprocess.DEVNULL)

    def _call(self, query_fp, database_fp, output_fp, **kwargs):
        """Call the BLAST program."""
        args = [
            "blastn",
            "-evalue", "1e-5",
            "-outfmt", "6 " + BLAST_FMT,
            ]
        for arg, val in kwargs.items():
            arg = "-" + arg
            if val is None:
                args.append(arg)
            else:
                args += [arg, str(val)]
        args += [
            "-query", query_fp,
            "-db", database_fp,
            "-out", output_fp,
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
                (hit['qseqid'], hit['qseq']),
                (hit['sseqid'], hit['sseq']))

        # We are going to need some repair or realignment.
        qseq = self.query_seqs[hit['qseqid']]
        assert(len(qseq) == hit['qlen'])
        sseq = self.ref_seqs[hit['sseqid']]
        assert(len(sseq) == hit['slen'])

        if self._needs_realignment(hit):
            aligned_qseq, aligned_sseq = align_semiglobal(qseq, sseq)
            return AlignedPair(
                (hit['qseqid'], aligned_qseq),
                (hit['sseqid'], aligned_sseq))

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
    def _needs_realignment(hit):
        more_to_the_left = (hit['qstart'] > 1) and \
                           (hit['sstart'] > 1)
        more_to_the_right = (hit['qend'] < hit['qlen']) and \
                            (hit['send'] < hit['slen'])
        return (more_to_the_left or more_to_the_right)

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

def align_semiglobal(qseq, sseq):
    alignment = pairwise2.align.globalms(
        sseq, qseq,
        5, -4, -10, -0.5, #match, mismatch, gapopen, gapextend
        penalize_end_gaps=False, one_alignment_only=True)
    subj_seq = alignment[0][0]
    query_seq = alignment[0][1]
    return query_seq, subj_seq
