import optparse
import os

from unassign.algorithm import Unassigner
from unassign.search_blast import BlastAligner
from unassign.parse import parse_fasta

def main(argv=None):
    p = optparse.OptionParser()
    p.add_option("--query_fp", help=(
        "Query sequences filepath (FASTA format) [REQUIRED]"))
    p.add_option("--output_dir", help=(
        "Output directory [REQUIRED]"))
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

    if not os.path.exists(opts.output_dir):
        os.mkdir(opts.output_dir)

    def mkfp(filename):
        return os.path.join(opts.output_dir, filename)

    a = BlastAligner(opts.type_strain_fp, opts.reference_fp)
    a.num_cpus = opts.num_cpus
    a.species_input_fp = mkfp("unassigner_query.fasta")
    a.species_output_fp = mkfp("unassigner_query_blastn.txt")
    a.refseq_input_fp = mkfp("unassigner_strains.fasta")
    a.refseq_output_fp = mkfp("unassigner_strains_blastn.txt")

    u = Unassigner(a)
    results = u.unassign(query_seqs)

    with open(mkfp("unassigner_output.tsv"), "w") as f:
        f.write(
            "QueryID\tTypestrainID\tRefseqID\t"
            "RegionMatch\tRegionTotal\tGlobalMatch\tGlobalTotal\n")
        for res in u.unassign(query_seqs):
            f.write("\t".join(map(str, res)))
            f.write("\n")
