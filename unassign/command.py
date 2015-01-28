import itertools
import logging
import optparse


class Unassigner(object):
    def __init__(self, aligner):
        self.aligner = aligner

    def unassign(self, query_seqs):
        """Execute unassignment algorithm on a filepath of query seqs.
        """
        logging.info("Aligning query seqs to type strain seqs")
        species_hits = self.aligner.search_species(query_seqs)

        logging.info("Aligning type strain seqs to reference seqs")
        refseq_hits = self.aligner.search_refseqs(species_hits)

        logging.info("Evaluating alignment results")
        grouped_refseq_hits = itertools.groupby(
            refseq_hits, key=lambda x: x.query_id)
        r_hits_species_id, r_hits = next(grouped_refseq_hits)
        for s_hit in species_hits:
            query_id = s_hit.query_id
            species_id = s_hit.subject_id
            if species_id != r_hits_species_id:
                r_hits_species_id, r_hits = next(grouped_refseq_hits)
            if species_id != r_hits_species_id:
                raise RuntimeError(
                    "Missing reference hits for %s" % species_id)
            for res in self._evaluate_confidence(s_hit, r_hits):
                yield res

    def _evaluate_confidence(self, species_hit, refseq_hits):
        """Compute probability of species attribution.
        """
        query_id = species_hit.query_id
        species_id = species_hit.subject_id
        start = species_hit.start_pos
        end = species_hit.end_pos
        for r_hit in refseq_hits:
            refseq_id = r_hit.subject_id
            a, b = r_hit.count_matches(start, end)
            c, d = r_hit.count_matches()
            yield query_id, species_id, a, b, refseq_id, c, d


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
