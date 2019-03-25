import itertools
import logging
import math

import scipy
import scipy.special
import scipy.misc

from unassign.search_blast import BlastSearch, BlastRefiner

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


class UnassignAligner(object):
    def __init__(self, species_fp):
        self.species_fp = species_fp
        self.species_max_hits = 1
        self.species_input_fp = None
        self.species_output_fp = None

    def search_species(self, seqs):
        """Search species typestrains for match to query sequences."""
        b = BlastSearch(self.species_fp)
        r = BlastRefiner(seqs, self.species_fp)
        hits = b.search(
            seqs, self.species_max_hits,
            self.species_input_fp, self.species_output_fp)
        for hit in hits:
            yield r.refine_hit(hit)


class BasicAlgorithm(UnassignerAlgorithm):
    def __init__(self, aligner):
        super(BasicAlgorithm, self).__init__(aligner)
        self.prior_alpha = 0.5
        self.prior_beta = 0.5
        self.species_threshold = 0.975

    def unassign(self, query_seqs):
        logging.info("Aligning query seqs to type strain seqs")
        species_hits = self.aligner.search_species(query_seqs)
        for h in species_hits:
            res = self._get_probability(h)
            yield res

    def _get_probability(self, species_hit):
        region_matches, region_positions = species_hit.count_matches()
        region_mismatches = region_positions - region_matches

        alpha = region_mismatches + self.prior_alpha
        beta = region_matches + self.prior_beta

        total_positions = species_hit.subject_len
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
        for line in super(BasicAlgorithm, self).format(results):
            yield line
