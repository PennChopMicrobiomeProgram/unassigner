import argparse
import sys

from unassigner.parse import parse_fasta, write_fasta
from unassigner.trim import (
    CompleteMatcher, PartialMatcher,
    deambiguate, reverse_complement,
    )

def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument(
        "--input_file", type=argparse.FileType("r"), default=sys.stdin)
    p.add_argument(
        "--output_file", type=argparse.FileType("w"), default=sys.stdout)
    p.add_argument(
        "--unmatched_output_file", type=argparse.FileType("w"),
    )
    p.add_argument("--query", required=True)
    p.add_argument(
        "--max_mismatch", type=int, default=2,
        help="Maximum number of mismatches in complete match")
    p.add_argument(
        "--min_partial", type=int, default=5,
        help="Minimum length of partial sequence match")
    p.add_argument(
        "--reverse_complement_query", action="store_true",
        help="Reverse complement the query seq before search")

    args = p.parse_args(argv)

    queryset = deambiguate(args.query)
    if args.reverse_complement_query:
        queryset = [reverse_complement(q) for q in queryset]

    seqs = parse_fasta(args.input_file)

    cm = CompleteMatcher(queryset, args.max_mismatch)
    unmatched_cm = []
    for seq_id, seq in seqs:
        m = cm.find_match(seq)
        if m is None:
            unmatched_cm.append((seq_id, seq))
        else:
            write_fasta(args.output_file, [(seq_id, seq)])

    if args.min_partial > 0:
        pm = PartialMatcher(queryset, args.min_partial)
        unmatched_pm = []
        for seq_id, seq in unmatched_cm:
            m = pm.find_match(seq)
            if m is None:
                unmatched_pm.append((seq_id, seq))
            else:
                write_fasta(args.output_file, [(seq_id, seq)])
    else:
        unmatched_pm = unmatched_cm

    if args.unmatched_output_file is not None:
        write_fasta(args.unmatched_output_file, unmatched_pm)


