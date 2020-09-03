import collections
import itertools
import math
import operator

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
    def load_database(cls, f):
        for line in f:
            line = line.rstrip()
            toks = line.split("\t")
            typestrain_id = toks[0]
            ref_seq_id = toks[1]
            mismatch_positions = [int(x) for x in toks[2:]]
            typestrain_mm = cls.db.get(typestrain_id, list)
            cls.db[typestrain_id].append(mismatch_positions)

    def __init__(self, alignment):
        self.alignment = alignment

    def unassign_threshold(self, min_id=0.975):
        pass

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

