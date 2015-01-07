#!/usr/bin/env python
from __future__ import division
import optparse
from cStringIO import StringIO

from unassign.search_blast import blast_to, hit_region_identity


def parse_fasta(f):
    f = iter(f)
    desc = f.next().strip()[1:]
    seq = StringIO()
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            yield desc, seq.getvalue()
            desc = line[1:]
            seq = StringIO()
        else:
            seq.write(line)
    yield desc, seq.getvalue()


def main(argv=None):
    p = optparse.OptionParser()
    p.add_option("--query_fp", help=(
        "Query sequences filepath (FASTA format) [REQUIRED]"))
    p.add_option("--type_strain_fp", default="species.fasta", help=(
        "Type strain sequences filepath (FASTA format + BLAST database) "
        "[default: %default]"))
    p.add_option("--reference_fp", default="gg_13_5.fasta", help=(
        "Reference sequence filepath (FASTA format + BLAST database) "
        "[default: %default]"))
    opts, args = p.parse_args(argv)
    
    f = open(opts.type_strain_fp)
    recs = parse_fasta(f)
    strain_seqs = dict((desc.split()[0], seq) for desc, seq in recs)
    f.close()

    query_blastout_fp = "unassigner_query_blastn.txt"
    query_hits = blast_to(opts.query_fp, opts.type_strain_fp, query_blastout_fp)
    top_query_hit = query_hits[0]

    top_query_hit_fp = "unassigner_top_hit.fasta"
    top_strain_id = top_query_hit['sseqid']
    with open(top_query_hit_fp, "w") as f:
        f.write(">%s\n%s\n" % (top_strain_id, strain_seqs[top_strain_id]))
    strain_blastout_fp = "unassigner_strain_blastn.txt"
    strain_hits = blast_to(top_query_hit_fp, opts.reference_fp, strain_blastout_fp, max_hits=100)

    start = top_query_hit['sstart']
    end = top_query_hit['send']
    for hit in strain_hits:
        a,b,c,d = hit_region_identity(hit, start, end)
        print top_strain_id, hit['sseqid'], a, b, c, d
    
if __name__ == "__main__":
    main()
