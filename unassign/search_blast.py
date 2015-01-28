from __future__ import division
import itertools
import subprocess
import tempfile

from unassign.parse import write_fasta

BLAST_FMT = (
    "qseqid sseqid pident length mismatch gapopen "
    "qstart qend sstart send qlen slen qseq sseq")
BLAST_FIELDS = BLAST_FMT.split()
BLAST_FIELD_TYPES = [
    str, str, float, int, int, int,
    int, int, int, int, int, int, str, str]


class BlastAlignment(object):
    def __init__(self, hit):
        self._hit = hit
        self.subject_id = hit['sseqid']
        self.query_id = hit['qseqid']
        self.start_pos = hit['sstart']
        self.end_pos = hit['send']

    def count_matches(self, start=None, end=None):
        """See docstring for hit_identity."""
        return hit_identity(self._hit, start, end)


class BlastAligner(object):
    """Align sequences with BLAST."""
    alignment_cls = BlastAlignment

    def __init__(self, num_threads=1):
        self.num_threads = num_threads

    def search_species(self, query_fp, species_fp, output_fp=None, max_hits=1):
        """Search species typestrains for match to query sequences."""
        if output_fp is None:
            outfile = tempfile.NamedTemporaryFile()
            output_fp = outfile.name
        self._call(
            query_fp, species_fp, output_fp,
            max_target_seqs=max_hits,
            num_threads=self.num_threads)
        return self._load(output_fp)

    def search_refseqs(self, query_seqs, refseqs_fp, input_fp=None,
                       output_fp=None, max_hits=100):
        """Search reference seqs for matches to species typestrains."""
        if input_fp is None:
            infile = tempfile.NamedTemporaryFile()
            input_fp = infile.name
        with open(input_fp, "w") as f:
            write_fasta(f, query_seqs)

        if output_fp is None:
            outfile = tempfile.NamedTemporaryFile()
            output_fp = outfile.name
        self._call(
            input_fp, refseqs_fp, output_fp,
            max_target_seqs=max_hits,
            num_threads=self.num_threads)
        return self._load(output_fp)

    @classmethod
    def _load(self, output_fp):
        """Load hits from an output file."""
        with open(output_fp) as f:
            hits = [self.alignment_cls(x) for x in self._parse(f)]
        return hits

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

    @staticmethod
    def _index(fasta_fp):
        return subprocess.check_call([
            "makeblastdb",
            "-dbtype", "nucl",
            "-in", fasta_fp,
            ])

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


def hit_identity(hit, start=None, end=None):
    """Count regional and total matches in BLAST hit.

    Parameters
    ----------
    hit : a dictionary representing the BLAST hit, must have the
        following keys: qseq, sseq, qstart, qend, qlen
    start : start position in query sequence
    end : end position in query sequence

    Returns
    -------
    tuple containing two ints
        Number of matching positions and total number of query
        nucleotides in the alignment.  Portions of the query sequence
        occurring over terminal gaps are not counted in the total.

    Notes
    -----
    Because sequence positions are indexed from 1 in BLAST, the
    start and end positions are indexed starting from 1.
    """
    if start is None:
        start = 1
    if end is None:
        end = hit['qlen']

    total = 0
    matches = 0

    # Count matches in alignment region.
    qpos = hit['qstart']
    for qchar, hchar in zip(hit['qseq'], hit['sseq']):
        if start <= qpos <= end:
            total += 1
            if qchar == hchar:
                matches += 1
        if qchar != '-':
            qpos += 1

    # If query is not aligned at the start position, count positions
    # from the start position to the left of the alignment as
    # mismatches.
    if start < hit['qstart']:
        # Number of nts in query from left end of alignment to query
        # start position.
        query_left = hit['qstart'] - start

        # Want to make sure there are at least this many nts available
        # for alignment in the subject sequence.  Number of nts in
        # subject to left of alignment.
        subject_left = hit['sstart'] - 1

        # Add the minimum of query_left vs. subject_left to the total.
        if query_left < subject_left:
            total += query_left
        else:
            total += subject_left

    # If query is not aligned at the end position, count positions to the
    # right of the alignment as mismatches.
    if end > hit['qend']:
        # Number of nts in query from right end of alignment to query
        # end position.
        query_right = end - hit['qend']

        # Number of nts in subject to right of alignment.
        subject_right = hit['slen'] - hit['send']

        # Again, add the smaller number to the total.
        if query_right < subject_right:
            total += query_right
        else:
            total += subject_right

        
    return matches, total

