import argparse
import os

from unassign.algorithm import (
    UnassignAligner, FileAligner, ThresholdAlgorithm,
)
from unassign.parse import parse_fasta

def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("query_fp", help=(
        "Query sequences filepath (FASTA format)"))
    p.add_argument("output_dir", help=(
        "Output directory"))
    p.add_argument("-t", "--type_strain_fp", default="species.fasta", help=(
        "Type strain sequences filepath (FASTA format + BLAST database) "
        "[default: %(default)s]"))
    p.add_argument("--alignment_fp", type=argparse.FileType("r"),
        help="Use pre-computed alignment file")
    p.add_argument("--verbose", action="store_true", help=(
        "Activate verbose mode."))
    args = p.parse_args(argv)

    if args.verbose is True:
        logging.basicConfig(level=logging.INFO)

    with open(args.query_fp) as f:
        query_seqs = list(parse_fasta(f, trim_desc=True))

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    def mkfp(filename):
        return os.path.join(args.output_dir, filename)

    if args.alignment_fp:
        a = FileAligner(args.type_strain_fp, args.alignment_fp)
    else:
        a = UnassignAligner(args.type_strain_fp)
        a.species_input_fp = mkfp("unassigner_query.fasta")
        a.species_output_fp = mkfp("unassigner_query_blastn.txt")

    u = ThresholdAlgorithm(a)
    results = u.unassign(query_seqs)

    with open(mkfp("unassigner_output.tsv"), "w") as f:
        for line in u.format(results):
            f.write(line)

