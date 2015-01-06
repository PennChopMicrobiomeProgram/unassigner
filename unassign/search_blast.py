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


def hit_region_identity(hit, qstart=0, qend=0):
    in_total = 0
    in_matches = 0

    total = 0
    matches = 0
    
    # matches in alignment region
    qpos = hit['qstart']
    for qchar, hchar in zip(hit['qseq'], hit['sseq']):
        if qstart <= qpos <= qend:
            in_total += 1
            if qchar == hchar:
                in_matches += 1
        total += 1
        if qchar == hchar:
            matches += 1
        if qchar != '-':
            qpos += 1

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
        
    return in_matches, in_total, matches, total

