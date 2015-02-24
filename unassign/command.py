import optparse

from unassign.algorithm import Unassigner
from unassign.search_blast import BlastAligner

def main(argv=None):
    p = optparse.OptionParser()
    p.add_option("--query_fp", help=(
        "Query sequences filepath (FASTA format) [REQUIRED]"))
    p.add_option("--type_strain_fp", default="species.fasta", help=(
        "Type strain sequences filepath (FASTA format + BLAST database) "
        "[default: %default]"))
    p.add_option("--reference_fp", default="refseqs.fasta", help=(
        "Reference sequence filepath (FASTA format + BLAST database) "
        "[default: %default]"))
    p.add_option("--num_cpus", type="int", default=1, help=(
        "Number of CPUs to use in seqrch and alignment steps "
        "[default: %default]"))
    p.add_option("--verbose", action="store_true", help=(
        "Activate verbose mode."))
    opts, args = p.parse_args(argv)

    if opts.verbose is True:
        logging.basicConfig(level=logging.INFO)

    with open(opts.query_fp) as f:
        query_seqs = list(parse_fasta(f, trim_desc=True))

    a = BlastAligner(opts.type_strain_fp, opts.reference_fp)
    a.num_cpus = opts.num_cpus
    a.species_input_fp = "unassigner_query.fasta"
    a.species_output_fp = "unassigner_query_blastn.txt"
    a.refseq_input_fp = "unassigner_strains.fasta"
    a.refseq_output_fp = "unassigner_strains_blastn.txt"

    u = Unassigner(a)
    for res in u.unassign(query_seqs):
        print "\t".join(res)
