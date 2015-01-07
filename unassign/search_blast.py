from __future__ import division
import subprocess

BLAST_FMT = (
    "qseqid sseqid pident length mismatch gapopen "
    "qstart qend sstart send qlen slen qseq sseq")
BLAST_FIELDS = BLAST_FMT.split()
BLAST_FIELD_TYPES = [
    str, str, float, int, int, int,
    int, int, int, int, int, int, str, str]

def check_call_blastn(query_fp, database_fp, output_fp, max_hits=100):
    args = [
        "blastn",
        "-evalue", "1e-5",
        "-outfmt", "6 " + BLAST_FMT,
        "-max_target_seqs", str(max_hits),
        "-query", query_fp,
        "-db", database_fp,
        "-out", output_fp,
        ]
    subprocess.check_call(args)


def parse_blastn(f):
    for line in f:
        line = line.strip()
        if line.startswith("#"):
            continue
        vals = line.split("\t")
        vals = [fn(v) for fn, v in zip(BLAST_FIELD_TYPES, vals)]
        yield dict(zip(BLAST_FIELDS, vals))


def blast_to(query_fp, database_fp, output_fp, max_hits=100):
    check_call_blastn(query_fp, database_fp, output_fp, max_hits)
    with open(output_fp) as f:
        return list(parse_blastn(f))


def hit_identity(hit):
    total = 0
    matches = 0

    # matches in alignment region
    for qchar, hchar in zip(hit['qseq'], hit['sseq']):
        total = total + 1
        if qchar == hchar:
            matches = matches + 1

    # if query is not aligned at first base
    if hit['qstart'] > 1:
        # Count unaligned nts as part of total
        total = total + hit['qstart'] - 1
        
        # TODO: check if we are at the end of subject sequence and do
        # not count if subject can't be extended.

    # if query is not aligned at last base
    if hit['qend'] < hit['qlen']:
        # Count unaligned nts as part of total
        total = total + hit['qlen'] - hit['qend']

        # TODO: again, check if we are out of range in the subject seqs
        
    return matches / total


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
        Number of matching positions and number of nucleotides in query.

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

