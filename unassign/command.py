import argparse
import os

from unassign.algorithm import (
    NoRefseqsAlgorithm, RefseqsAlgorithm,
    )
from unassign.search_blast import BlastAligner
from unassign.parse import parse_fasta

def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("query_fp", help=(
        "Query sequences filepath (FASTA format)"))
    p.add_argument("output_dir", help=(
        "Output directory"))
    p.add_argument("--noref", action="store_true", help=(
        "Use simple estimate with no reference sequences"))
    p.add_argument("--type_strain_fp", default="species.fasta", help=(
        "Type strain sequences filepath (FASTA format + BLAST database) "
        "[default: %(default)s]"))
    p.add_argument("--reference_fp", default="refseqs.fasta", help=(
        "Reference sequence filepath (FASTA format + BLAST database) "
        "[default: %(default)s]"))
    p.add_argument("--num_cpus", type=int, default=1, help=(
        "Number of CPUs to use in seqrch and alignment steps "
        "[default: %(default)s]"))
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

    a = BlastAligner(args.type_strain_fp, args.reference_fp)
    a.num_cpus = args.num_cpus
    a.species_input_fp = mkfp("unassigner_query.fasta")
    a.species_output_fp = mkfp("unassigner_query_blastn.txt")
    a.refseq_input_fp = mkfp("unassigner_strains.fasta")
    a.refseq_output_fp = mkfp("unassigner_strains_blastn.txt")

    if args.noref:
        u = NoRefseqsAlgorithm(a)
    else:
        u = RefseqsAlgorithm(a)
    results = u.unassign(query_seqs)

    with open(mkfp("unassigner_output.tsv"), "w") as f:
        for line in u.format(results):
            f.write(line)

