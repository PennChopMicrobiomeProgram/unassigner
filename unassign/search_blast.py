from __future__ import division
import itertools
import subprocess
import tempfile

from unassign.parse import write_fasta, load_fasta
from unassign.util import uniq, count_while_equal
from unassign.alignment import Alignment, Aligner



class BlastAlignment(Alignment):
    def __init__(self, hit):
        self._hit = hit
        self.subject_id = hit['sseqid']
        self.query_id = hit['qseqid']
        self.query_seq, self.subject_seq = self._get_aligned_seqs(hit)
        assert(len(self.query_seq) == len(self.subject_seq))

    def alignment_length(self):
        return _hit_identity(self._hit)[1]

    def num_matches(self):
        return _hit_identity(self._hit)[0]

    def unaligned_length(self):
        total_length = len(
            [x for x in self.subject_seq if x != '-'])
        return total_length - self.alignment_length()

    @property
    def start_idx(self):
        return max(
            count_while_equal(self.query_seq, "-"),
            count_while_equal(self.subject_seq, "-"))

    @property
    def end_idx(self):
        return len(self.query_seq) - max(
            count_while_equal(reversed(self.query_seq), "-"),
            count_while_equal(reversed(self.subject_seq), "-"))

    def get_pairs(self, start=None, end=None):
        astart = self.start_idx
        if (start is None) or (start < astart):
            start = astart

        aend = self.end_idx
        if (end is None) or (end > aend):
            end = aend

        return list(zip(self.query_seq, self.subject_seq))[start:end]

    def _get_aligned_seqs(self, hit):
        # Number of nts outside the local alignment
        query_nleft = hit['qstart'] - 1
        subj_nleft = hit['sstart'] - 1
        query_nright = hit['qlen'] - hit['qend']
        subj_nright = hit['slen'] - hit['send']

        # Number of positions to left and right of local alignment.
        nleft = max(query_nleft, subj_nleft)
        nright = max(query_nright, subj_nright)

        # Fill in query outside the local alignment with "X"
        query_lgaps = "-" * (nleft - query_nleft)
        query_lfill = "X" * query_nleft
        query_rfill = "X" * query_nright
        query_rgaps = "-" * (nright - query_nright)
        query_seq = ''.join(
            query_lgaps + query_lfill +
            hit['qseq'] +
            query_rfill + query_rgaps)

        # Fill in subj outside the local alignment with "H"
        subj_lgaps = "-" * (nleft - subj_nleft)
        subj_lfill = "H" * subj_nleft
        subj_rfill = "H" * subj_nright
        subj_rgaps = "-" * (nright - subj_nright)
        subj_seq = ''.join(
            subj_lgaps + subj_lfill +
            hit['sseq'] +
            subj_rfill + subj_rgaps)

        return query_seq, subj_seq

    def get_local_pairs(self):
        return zip(self._hit['qseq'], self._hit['sseq'])
    
    def count_matches(self, start=None, end=None):
        """See docstring for _hit_identity."""
        return _hit_identity(self._hit, start, end)


class BlastAligner(Aligner):
    """Align sequences with BLAST."""
    alignment_cls = BlastAlignment
    executable = "blastn"
    blastfmt = (
        "qseqid sseqid pident length mismatch gapopen "
        "qstart qend sstart send qlen slen qseq sseq")
    field_names = blastfmt.split()
    field_types = [
        str, str, float, int, int, int,
        int, int, int, int, int, int, str, str]

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
            "-outfmt", "6 " + self.blastfmt,
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


def _hit_identity(hit, start=None, end=None):
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

