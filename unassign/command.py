import logging
import optparse
import sys


from unassign.search_blast import (
    blast_to, top_hits, hit_identity, group_by_query,
    )
from unassign.parse import parse_fasta, write_fasta, load_fasta
from unassign.util import uniq


class Unassigner(object):
    def __init__(self, species_fp, refseqs_fp, aligner):
        self.species_fp = species_fp
        self.refseqs_fp = refseqs_fp
        self.aligner = aligner

    def unassign(self, query_fp):
        """Execute unassignment algorithm on a filepath of query seqs.
        """
        logging.info("Aligning query seqs to type strain seqs")
        species_hits = self.aligner.search_species(
            query_fp, self.species_fp,
            output_fp="unassigner_query_blastn.txt")

        logging.info("Aligning type strain seqs to reference seqs")
        species_seqs = self._get_species_seqs(species_hits)
        refseq_hits = self.aligner.search_refseqs(
            species_seqs, self.refseqs_fp,
            input_fp="unassigner_top_hit.fasta",
            output_fp="unassigner_strain_blastn.txt")

        logging.info("Evaluating aignment results")
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
            self._evaluate_confidence(s_hit, r_hits)

    def _get_species_seqs(self, hits):
        """Fetch seqs for each species in the list of hits.

        Each unique species in the list is returned only once in the
        list of results, but the order of hits is preserved to
        facilitate caching.
        """
        species_ids = uniq(x.subject_id for x in species_hits)
        seqs = load_fasta(self.species_fp)
        return [(x, seqs[x]) for x in species_ids]

    def _evaluate_confidence(self, species_hit, refseq_hits):
        """Compute probability of species attribution.
        """
        query_id = species_hit['qseqid']
        species_id = species_hit['sseqid']
        start = species_hit['sstart']
        end = species_hit['send']
        for r_hit in refseq_hits:
            refseq_id = r_hit['sseqid']
            a, b = hit_identity(r_hit, start, end)
            c, d = hit_identity(r_hit)
            print query_id, species_id, a, b, refseq_id, c, d


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

    a = BlastAligner(opts.num_cpus)
    u = Unassigner(opts.type_strain_fp, opts.reference_fp, a)
    u.unassign(opts.query_fp)
