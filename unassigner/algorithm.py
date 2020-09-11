import collections
import itertools
import math
import operator

import numpy
import scipy
import scipy.special

from unassigner.alignment import AlignedRegion
from unassigner.align import VsearchAligner, HitExtender
from unassigner.parse import parse_fasta

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

class VariableMismatchRate:
    """Predict unobserved mismatches by estimating a mismatch rate
    """
    db = collections.defaultdict(list)

    @classmethod
    def clear_database(cls):
        cls.db.clear()

    @classmethod
    def load_database(cls, f):
        for line in f:
            line = line.rstrip()
            toks = line.split("\t")
            typestrain_id = toks[0]
            ref_seq_id = toks[1]
            mismatch_positions = [int(x) for x in toks[2:]]
            cls.db[typestrain_id].append(mismatch_positions)

    @classmethod
    def _get_mismatches(self, typestrain_id, start_idx, end_idx):
        refs = self.db[typestrain_id]
        for ref_positions in refs:
            mismatch_is_in_region = [
                (pos >= start_idx) & (pos <= end_idx)
                for pos in ref_positions]
            region_mms = mismatch_is_in_region.count(True)
            nonregion_mms = mismatch_is_in_region.count(False)
            yield (region_mms, nonregion_mms)

    def __init__(self, alignment):
        self.alignment = alignment
        self.query_id = alignment.query_id

    def unassign_threshold(self, min_id=0.975):
        # Use all the beta-binomial logic from ConstantMismatchRate,
        # just adjust alpha and beta based on reference
        # sequences. Here's how. Reparameterize beta as mu and v,
        # following wikipedia. Hold v constant. We are going to update
        # mu. In the constant rate algorithm, mu2 = mu1. In the
        # variable rate algorithm, we determine log(mu2 / mu1) by
        # averaging the observed values from the reference
        # sequences. To stabilize things, start with a list of
        # [0,0,0,0,0]. Then, for each reference sequence, compute mu2
        # and mu1, take the log, and append to the list. Average the
        # values in the list. Now use this as the new value of mu2 for
        # the query sequence.

        # Clip out the aligned region
        region = AlignedRegion.without_endgaps(self.alignment)
        region_alignment = region.trim_ends()
        region_positions = region_alignment.alignment_len
        region_matches = region_alignment.count_matches()
        region_mismatches = region_positions - region_matches
        region_subject_positions = region_alignment.subject_len

        # Calcuate alpha, beta, mu, and v in aligned region
        alpha1 = region_mismatches + 0.5
        beta1 = region_matches + 0.5
        v1 = alpha1 + beta1
        mu1 = alpha1 / v1

        # Compute number of positions outside aligned region
        nonregion_subject_positions = (
            self.alignment.subject_len - region_subject_positions)
        total_positions = (
            region_positions + nonregion_subject_positions)

        # Get mismatches from database
        typestrain_id = self.alignment.subject_id
        typestrain_start_idx, typestrain_end_idx = region.in_subject()
        reference_mismatches = self._get_mismatches(
            typestrain_id, typestrain_start_idx, typestrain_end_idx)

        # Get estimate for gamma = log(mu2 / mu1)
        # From reference alignments
        reference_logvals = [0, 0, 0, 0, 0]
        for region_mms, nonregion_mms in reference_mismatches:
            # mu = alpha / (alpha + beta)
            # alpha = mismatches + 0.5
            # beta = matches + 0.5
            # matches = len - mismatches
            # beta = len - mismatches + 0.5
            # mu = (mismatches + 0.5) / (len + 1)
            ref_mu1 = (region_mms + 0.5) / (region_subject_positions + 1)
            ref_mu2 = (nonregion_mms + 0.5) / (nonregion_subject_positions + 1)
            log_mu2_mu1 = math.log(ref_mu2 / ref_mu1)
            reference_logvals.append(log_mu2_mu1)
        # TODO: add weighting
        gamma = numpy.mean(reference_logvals)

        # Calculate mu2, get alpha2 and beta2
        # log(mu2 / mu1) = gamma
        # log(mu2) - log(mu1) = gamma
        # log(mu2) = log(mu1) + gamma
        # mu2 = exp(log(mu1) + gamma)
        mu2 = math.exp(math.log(mu1) + gamma)
        v2 = v1
        # alpha2, beta2
        # mu2 = alpha2 / v2
        # alpha2 = mu2 * v2
        alpha2 = mu2 * v2
        # v2 = alpha2 + beta2
        beta2 = v2 - alpha2

        # Maximum number of mismatches outside observed region
        species_mismatch_threshold = 1 - min_id
        max_total_mismatches = int(math.floor(
            species_mismatch_threshold * total_positions))
        max_nonregion_mismatches = max_total_mismatches - region_mismatches

        # Compute probability
        prob_compatible = beta_binomial_cdf(
            max_nonregion_mismatches, nonregion_subject_positions,
            alpha2, beta2)
        prob_incompatible = 1 - prob_compatible

        return {
            "typestrain_id": self.alignment.subject_id,
            "region_mismatches": region_mismatches,
            "region_positions": region_positions,
            "probability_incompatible": prob_incompatible,
            "mu1": mu1,
            "num_references": len(reference_logvals),
            "mu2": mu2,
            "nonregion_positions_in_subject": nonregion_subject_positions,
            "max_nonregion_mismatches": max_nonregion_mismatches,
        }


class ConstantMismatchRate:
    """Predict unobserved mismatches assuming a constant mismatch rate
    """

    result_keys = [
        "typestrain_id", "probability_incompatible", "region_mismatches",
        "region_positions", "region_matches", "nonregion_positions_in_subject",
        "max_nonregion_mismatches",
    ]
    null_result = dict((key, "NA") for key in result_keys)
    prior_alpha = 0.5
    prior_beta = 0.5

    def __init__(self, alignment):
        self.alignment = alignment
        self.query_id = alignment.query_id
        self.region = alignment.trim_endgaps()
        self.region_positions = self.region.alignment_len
        self.region_matches = self.region.count_matches()
        self.region_mismatches = self.region_positions - self.region_matches

        self.alpha = self.region_mismatches + self.prior_alpha
        self.beta = self.region_matches + self.prior_beta

    def unassign_threshold(self, min_id=0.975):
        """Unassign with a hard threshold

        Here, we set a threshold value for sequence similarity to the
        type strain sequence.  For a query sequence, we calculate the
        probability of falling below this similarity threshold over
        the full length of the 16S gene.  This value is the
        unassignment probability.
        """
        nonregion_subject_positions = (
            self.alignment.subject_len - self.region.subject_len)
        total_positions = (
            self.region_positions + nonregion_subject_positions)

        species_mismatch_threshold = 1 - min_id
        max_total_mismatches = int(math.floor(
            species_mismatch_threshold * total_positions))
        max_nonregion_mismatches = max_total_mismatches - self.region_mismatches

        prob_compatible = beta_binomial_cdf(
            max_nonregion_mismatches, nonregion_subject_positions,
            self.alpha, self.beta)
        prob_incompatible = 1 - prob_compatible

        return {
            "typestrain_id": self.alignment.subject_id,
            "probability_incompatible": prob_incompatible,
            "region_mismatches": self.region_mismatches,
            "region_positions": self.region_positions,
            "region_matches": self.region_matches,
            "nonregion_positions_in_subject": nonregion_subject_positions,
            "max_nonregion_mismatches": max_nonregion_mismatches,
        }


class UnassignerApp:
    def __init__(self, aligner, mm_rate = ConstantMismatchRate):
        self.aligner = aligner
        self.mm_rate = mm_rate
        self.alignment_min_percent_id = 0.975

    def unassign(self, query_seqs):
        query_seqs = list(query_seqs)
        query_ids = [seq_id for seq_id, seq in query_seqs]

        # Step 1.
        # Align query sequences to type strain sequences. Same for all
        # algorithms.
        alignments = self._align_query_to_type_strain(query_seqs)

        # Step 2.
        # For each query-type strain alignment, estimate distribution of
        # mismatches outside fragment. Different for constant vs. variable
        # mismatch algorithms.
        mm_rates = [self.mm_rate(a) for a in alignments]

        # Step 3.
        # For each query-type strain alignment, estimate unassignment
        # probability. Different for hard vs. soft threshold algorithms.
        results = [(r.query_id, r.unassign_threshold()) for r in mm_rates]

        # Step 4.
        # Group by query and yield results to caller. Same for all algorithms.
        results_by_query = collections.defaultdict(list)
        for query_id, res in results:
            results_by_query[query_id].append(res)
        for query_id in query_ids:
            query_results = results_by_query[query_id]
            if not query_results:
                query_results = [self.mm_rate.null_result]
            yield query_id, query_results

    def _align_query_to_type_strain(self, query_seqs):
        # We expect query_seqs to be a list
        query_ids = [seq_id for seq_id, seq in query_seqs]
        unsorted_alignments = self.aligner.search_species(query_seqs)

        alignments_by_query = collections.defaultdict(list)
        for a in unsorted_alignments:
            alignments_by_query[a.query_id].append(a)

        for query_id in query_ids:
            query_alignments = alignments_by_query[query_id]
            query_alignments = self._filter_alignments(query_alignments)
            for a in query_alignments:
                yield a

    def _filter_alignments(self, query_alignments):
        sorted_alignments = list(sorted(
            query_alignments, key=operator.attrgetter('percent_id'),
            reverse=True))
        filtered_alignments = [
            a for a in sorted_alignments
            if a.percent_id > self.alignment_min_percent_id]
        # Return one low-identity result if we have nothing better
        if sorted_alignments and not filtered_alignments:
            return sorted_alignments[:1]
        return filtered_alignments

