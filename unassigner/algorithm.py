import collections
import itertools
import math
import operator

import scipy
import scipy.special

from unassigner.alignment import AlignedRegion
from unassigner.align import VsearchAligner, HitExtender
from unassigner.parse import parse_fasta

class UnassignerAlgorithm:
    def __init__(self, aligner):
        self.aligner = aligner

    def unassign(self, query_seqs):
        query_seqs = list(query_seqs)
        query_ids = [seq_id for seq_id, seq in query_seqs]
        alignments_by_query_id = collections.defaultdict(list)

        for alignment in self.aligner.search_species(query_seqs):
            alignments_by_query_id[alignment.query_id].append(alignment)

        for query_id in query_ids:
            alignments = alignments_by_query_id[query_id]
            results = self._get_probability(alignments)
            yield query_id, list(results)


def beta_binomial_pdf(k, n, alpha, beta):
    binom_coeff = scipy.special.comb(n, k)
    if binom_coeff == 0:
        return 0
    t1 = math.log(binom_coeff)
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
        self.species_input_fp = None
        self.species_output_fp = None
        self.num_cpus = None

    def search_species(self, query_seqs):
        b = VsearchAligner(self.species_fp)
        vsearch_args = {
            "min_id": 0.9,
            "maxaccepts": 5,
        }
        if self.num_cpus:
            vsearch_args["threads"] = self.num_cpus
        hits = b.search(
            query_seqs,
            self.species_input_fp, self.species_output_fp, **vsearch_args)

        with open(self.species_fp) as f:
            ref_seqs = list(parse_fasta(f, trim_desc=True))
        xt = HitExtender(query_seqs, ref_seqs)
        for hit in hits:
            yield xt.extend_hit(hit)


class FileAligner:
    def __init__(self, species_fp, output_fp):
        self.species_fp = species_fp
        self.output_fp = output_fp

    def search_species(self, seqs):
        with open(self.species_fp) as f:
            ref_seqs = list(parse_fasta(f, trim_desc=True))
        xt = HitExtender(seqs, ref_seqs)
        with open(self.output_fp) as of:
            hits = VsearchAligner._parse(of)
            for hit in hits:
                yield xt.extend_hit(hit)


class ThresholdAlgorithm(UnassignerAlgorithm):
    """Threshold algorithm for species unassignment

    In this algorithm, we set a threshold value for sequence
    similarity to the type strain sequence.  For a query sequence, we
    calculate the probability of falling below this similarity
    threshold over the full length of the 16S gene.  This value is the
    unassignment probability.
    """
    result_keys = [
        "typestrain_id", "probability_incompatible", "region_mismatches",
        "region_positions", "region_matches", "nonregion_positions_in_subject",
        "max_nonregion_mismatches",
    ]
    null_result = dict((key, "NA") for key in result_keys)

    def __init__(self, aligner):
        super().__init__(aligner)
        self.prior_alpha = 0.5
        self.prior_beta = 0.5
        self.species_threshold = 0.975

    def _get_probability(self, hits):
        res = [self._get_indiv_probability(hit) for hit in hits]
        # Return something even if no hits were obtained
        if not res:
            return [self.null_result]
        res = list(
            sorted(res, key=operator.itemgetter('probability_incompatible')))
        res_filtered = [
            r for r in res if self.sequence_id(r) > self.species_threshold]
        # Return one low-identity result if we have nothing better
        if not res_filtered:
            return res[:1]
        # Otherwise return hits exceeding the identity threshold
        return res_filtered

    @staticmethod
    def sequence_id(result):
        return result["region_matches"] / result["region_positions"]

    def _get_indiv_probability(self, alignment):
        region = AlignedRegion.without_endgaps(alignment).trim_ends()
        region_positions = region.alignment_len
        region_matches = region.count_matches()
        region_mismatches = region_positions - region_matches

        alpha = region_mismatches + self.prior_alpha
        beta = region_matches + self.prior_beta

        nonregion_subject_positions = (
            alignment.subject_len - region.subject_len)
        total_positions = (
            region_positions + nonregion_subject_positions)

        species_mismatch_threshold = 1 - self.species_threshold
        max_total_mismatches = int(math.floor(
            species_mismatch_threshold * total_positions))
        max_nonregion_mismatches = max_total_mismatches - region_mismatches

        prob_compatible = beta_binomial_cdf(
            max_nonregion_mismatches, nonregion_subject_positions, alpha, beta)
        prob_incompatible = 1 - prob_compatible

        return {
            "typestrain_id": alignment.subject_id,
            "probability_incompatible": prob_incompatible,
            "region_mismatches": region_mismatches,
            "region_positions": region_positions,
            "region_matches": region_matches,
            "nonregion_positions_in_subject": nonregion_subject_positions,
            "max_nonregion_mismatches": max_nonregion_mismatches,
        }
