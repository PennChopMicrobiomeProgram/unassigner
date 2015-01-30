import itertools
import logging


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


