import itertools
import logging
import math

import scipy
import scipy.special
import scipy.misc


class UnassignerAlgorithm(object):
    def __init__(self, aligner):
        self.aligner = aligner

    def unassign(self, query_seqs):
        """Execute unassignment algorithm on query seqs."""
        raise NotImplemented

    def format(self, results):
        """Format results for output file."""
        for r in results:
            yield '\t'.join(map(str, r)) + "\n"


def beta_binomial_pdf(k, n, alpha, beta):
    t1 = math.log(scipy.misc.comb(n, k))
    t2 = scipy.special.betaln(k + alpha, n - k + beta)
    t3 = scipy.special.betaln(alpha, beta)
    logf = t1 + t2 - t3
    return math.exp(logf)


def beta_binomial_cdf(k_max, n, alpha, beta):
    k = 0
    val = 0
    while k <= k_max:
        val += beta_binomial_pdf(k, n, alpha, beta)
        k += 1
    return val


class NoRefseqsAlgorithm(UnassignerAlgorithm):
    def __init__(self, aligner):
        super(NoRefseqsAlgorithm, self).__init__(aligner)
        self.prior_alpha = 0.5
        self.prior_beta = 0.5
        self.species_threshold = 0.975

    def unassign(self, query_seqs):
        logging.info("Aligning query seqs to type strain seqs")
        species_hits = self.aligner.search_species(query_seqs)

        for h in species_hits:
            res = self._evaluate_noref_probability(h)
            yield res

    def _evaluate_noref_probability(self, species_hit):
        region_matches, region_positions = species_hit.count_matches()
        region_mismatches = region_positions - region_matches

        alpha = region_mismatches + self.prior_alpha
        beta = region_matches + self.prior_beta

        total_positions = len(
            [x for x in species_hit.subject_seq if x != '-'])
        nonregion_positions = total_positions - region_positions

        species_mismatch_threshold = 1 - self.species_threshold
        max_total_mismatches = int(math.floor(
            species_mismatch_threshold * total_positions))
        max_nonregion_mismatches = max_total_mismatches - region_mismatches

        prob_compatible = beta_binomial_cdf(
            max_nonregion_mismatches, nonregion_positions, alpha, beta)
        prob_new_species = 1 - prob_compatible

        return (species_hit.query_id, species_hit.subject_id,
            region_mismatches, region_matches,
            nonregion_positions, max_nonregion_mismatches,
            prob_new_species)

    def format(self, results):
        yield (
            "QueryID\tTypestrainID\t"
            "RegionMismatches\tRegionMatches\t"
            "NonregionPositions\tMaxNonregionMismatches\t"
            "ProbabilityNotThisSpecies\n")
        for line in super(NoRefseqsAlgorithm, self).format(results):
            yield line


class RefseqsAlgorithm(UnassignerAlgorithm):
    def unassign(self, query_seqs):
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
        start = species_hit.start_idx
        end = species_hit.end_idx
        for r_hit in refseq_hits:
            refseq_id = r_hit.subject_id
            a, b = r_hit.count_matches(start, end)
            c, d = r_hit.count_matches()
            yield query_id, species_id, a, b, refseq_id, c, d

    def format(self, results):
        yield (
            "QueryID\tTypestrainID\tRefseqID\t"
            "RegionMatch\tRegionTotal\t"
            "GlobalMatch\tGlobalTotal\n")
        for line in super(RefseqsAlgorithm, self).format(results):
            yield line

