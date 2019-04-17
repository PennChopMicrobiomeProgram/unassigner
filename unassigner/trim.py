import abc
import argparse
import collections
import itertools
import os
import shutil
import sys
import tempfile

from unassigner.parse import parse_fasta, write_fasta
from unassigner.align import VsearchAligner, HitExtender
from unassigner.alignment import AlignedRegion

class TrimmableSeqs:
    def __init__(self, recs):
        replicate_seqs = collections.defaultdict(list)
        self.descs = dict()
        for desc, seq in recs:
            seq_id = desc.split()[0]
            replicate_seqs[seq].append(seq_id)
            self.descs[seq_id] = desc

        self.seqs = dict()
        self.seq_ids = dict()
        for seq, seq_ids in replicate_seqs.items():
            rep_seq_id = seq_ids[0]
            self.seqs[rep_seq_id] = seq
            self.seq_ids[rep_seq_id] = seq_ids
        self.matches = dict()

    def all_matched(self):
        return set(self.matches) == set(self.seq_ids)

    def get_matched_offset0(self):
        for rep_seq_id, matchobj in self.matches.items():
            if matchobj.offset == 0:
                seq = self.seqs[rep_seq_id]
                yield rep_seq_id, seq

    def get_unmatched_recs(self):
        for rep_seq_id in self.seq_ids:
            if rep_seq_id not in self.matches:
                seq = self.seqs[rep_seq_id]
                yield rep_seq_id, seq

    def get_replicate_ids(self, rep_seq_id):
        return self.seq_ids[rep_seq_id]

    def get_desc(self, seq_id):
        return self.descs[seq_id]

    def get_replicate_recs(self, rep_seq_id):
        seq_ids = self.seq_ids[rep_seq_id]
        seq = self.seqs[rep_seq_id]
        for seq_id in seq_ids:
            yield seq_id, seq

    def register_match(self, rep_seq_id, matchobj):
        self.matches[rep_seq_id] = matchobj

    @classmethod
    def from_fasta(cls, f):
        recs = parse_fasta(f)
        return cls(recs)


PrimerMatch = collections.namedtuple(
    "PrimerMatch", ["start", "end", "offset", "message"])


class Matcher(abc.ABC):
    def __init__(self, queryset):
        self.queryset = queryset

    # TODO: parallelize this
    def find_in_seqs(self, seqs):
        recs = seqs.get_unmatched_recs()
        for seq_id, seq in recs:
            match = self.find_match(seq)
            if match is not None:
                seqs.register_match(seq_id, match)
                yield seq_id, match

    def find_match(self, seq):
        raise NotImplemented()


class CompleteMatcher(Matcher):
    def __init__(self, queryset, max_mismatch):
        super().__init__(queryset)
        self.max_mismatch = max_mismatch

        # To look for near matches in the reference sequence, we
        # generate the set of qeury sequences with 0, 1, or 2 errors
        # and look for these in the reference.  This is our
        # "mismatched queryset."
        possible_mismatches = range(max_mismatch + 1)
        self.mismatched_queryset = [
            self._mismatched_queries(n) for n in possible_mismatches]

    def _mismatched_queries(self, n_mismatch):
        # The generator is provides a sequence for one-time use, but
        # we need to go through it multiple times.  This function
        # wraps the generator to provide a list.
        return list(self._iter_mismatched_queries(n_mismatch))

    def _iter_mismatched_queries(self, n_mismatch):
        # This algorithm is terrible unless the number of mismatches is very small
        assert(n_mismatch in [0, 1, 2, 3])
        for query in self.queryset:
            idx_sets = itertools.combinations(range(len(query)), n_mismatch)
            for idx_set in idx_sets:
                # Replace base at each position with a literal "N", to match
                # ambiguous bases in the reference
                yield replace_with_n(query, idx_set)
                # Change to list because strings are immutable
                qchars = list(query)
                # Replace the base at each mismatch position with an
                # ambiguous base specifying all possibilities BUT the one
                # we see.
                for idx in idx_set:
                    qchars[idx] = AMBIGUOUS_BASES_COMPLEMENT[qchars[idx]]
                    # Expand to all possibilities for mismatching at this
                    # particular set of positions
                for query_with_mismatches in deambiguate(qchars):
                    yield query_with_mismatches

    def find_match(self, seq):
        for n_mismatches, queryset in enumerate(self.mismatched_queryset):
            for query in queryset:
                start_idx = seq.find(query)
                if start_idx > -1:
                    if n_mismatches == 0:
                        msg = "Exact"
                    elif n_mismatches == 1:
                        msg = "Complete, 1 mismatch"
                    else:
                        msg = "Complete, {0} mismatches".format(n_mismatches)
                    end_idx = start_idx + len(query)
                    obs_primer = seq[start_idx:end_idx]
                    return PrimerMatch(start_idx, end_idx, 0, msg)


def replace_with_n(seq, idxs):
    chars = list(seq)
    for idx in idxs:
        chars[idx] = "N"
    return "".join(chars)


class PartialMatcher(Matcher):
    def __init__(self, queryset, min_length):
        super().__init__(queryset)
        self.min_length = min_length

        self.partial_queries = set()
        for query in self.queryset:
            for partial_query in partial_seqs(query, self.min_length):
                self.partial_queries.add(partial_query)

    def find_match(self, seq):
        for partial_query in self.partial_queries:
            if seq.startswith(partial_query):
                end_idx = len(partial_query)
                return PrimerMatch(0, end_idx, 0, "Partial")


class AlignmentMatcher(Matcher):
    def __init__(
            self, alignment_dir, min_pct_id=75, min_aligned_frac=0.6,
            cores=1, suffix="0"):
        assert(os.path.exists(alignment_dir))
        assert(os.path.isdir(alignment_dir))
        self.alignment_dir = alignment_dir
        self.min_pct_id = min_pct_id
        self.min_aligned_frac = min_aligned_frac
        self.cores = cores
        self.suffix = suffix

    def _make_fp(self, filename):
        return os.path.join(self.alignment_dir, filename)

    def find_in_seqs(self, seqs):
        if seqs.all_matched():
            return

        # Create the file paths
        subject_fp = self._make_fp("subject_{0}.fa".format(self.suffix))
        query_fp = self._make_fp("query_{0}.fa".format(self.suffix))
        result_fp = self._make_fp("query_{0}.txt".format(self.suffix))

        # Search
        with open(subject_fp, "w") as f:
            write_fasta(f, seqs.get_matched_offset0())
        ba = VsearchAligner(subject_fp)
        search_args = {
            "min_id": round(self.min_pct_id / 100, 2),
            "top_hits_only": None}
        if self.cores > 0:
            search_args["threads"] = self.cores
        hits = ba.search(
            seqs.get_unmatched_recs(), input_fp=query_fp, output_fp=result_fp,
            **search_args)

        # Refine
        bext = HitExtender(seqs.get_unmatched_recs(), seqs.get_matched_offset0())
        for hit in hits:
            alignment = bext.extend_hit(hit)
            subject_match = seqs.matches[alignment.subject_id]
            aligned_region = AlignedRegion.from_subject(
                alignment, subject_match.start, subject_match.end)
            query_start_idx, query_end_idx = aligned_region.in_query()
            query_offset = aligned_region.query_offset()
            matchobj = PrimerMatch(
                query_start_idx, query_end_idx, query_offset, "Alignment")
            yield alignment.query_id, matchobj


def partial_seqs(seq, min_length):
    max_start_idx = len(seq) - min_length + 1
    for start_idx in range(1, max_start_idx):
        yield seq[start_idx:]

def aligned_frac(hit):
    unaligned_left = min(hit["qstart"], hit["sstart"])
    unaligned_right = min(
        hit["qlen"] - hit["qend"],
        hit["slen"] - hit["send"],
    )
    unaligned = unaligned_left + unaligned_right
    aligned = hit["length"]
    return aligned / (aligned + unaligned)

class TrimraggedApp(object):
    def __init__(self, seqs, trim_right, writer, min_trimmed_length):
        self.seqs = seqs
        self.trim_right = trim_right
        self.writer = writer
        self.min_trimmed_length = min_trimmed_length
        self.matchers = []

    def run(self):
        for matcher in self.matchers:
            matches = matcher.find_in_seqs(self.seqs)
            for rep_seq_id, matchobj in matches:
                if matchobj is not None:
                    self.seqs.register_match(rep_seq_id, matchobj)
                    seq_ids = self.seqs.seq_ids[rep_seq_id]
                    seq = self.seqs.seqs[rep_seq_id]
                    if self.trim_right:
                        trimmed_seq = trim_right(seq, matchobj)
                    else:
                        trimmed_seq = trim_left(seq, matchobj)
                    for seq_id in seq_ids:
                        desc = self.seqs.get_desc(seq_id)
                        if len(trimmed_seq) >= self.min_trimmed_length:
                            self.writer.write_trimmed(desc, trimmed_seq)
                        self.writer.write_stats(seq_id, seq, matchobj)

        for rep_seq_id, seq in self.seqs.get_unmatched_recs():
            for seq_id in self.seqs.get_replicate_ids(rep_seq_id):
                desc = self.seqs.get_desc(seq_id)
                self.writer.write_stats(seq_id, seq, None)


class Writer(object):
    def __init__(self, trimmed_file, stats_file):
        self.trimmed_file = trimmed_file
        self.stats_file = stats_file

    def write_trimmed(self, desc, seq):
        self.trimmed_file.write(">{0}\n{1}\n".format(desc, seq))

    def write_stats(self, seq_id, seq, matchobj):
        if matchobj is None:
            self.stats_file.write("{0}\tUnmatched\tNA\tNA\tNA\t\n".format(seq_id))
        else:
            matched_seq = seq[matchobj.start:matchobj.end]
            self.stats_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                seq_id, matchobj.message, matchobj.start, matchobj.end,
                matchobj.offset, matched_seq))


def trim_left(seq, matchobj):
    return seq[matchobj.end:]

def trim_right(seq, matchobj):
    return seq[:matchobj.start]

def trim_middle(seq, matchobj):
    return seq[matchobj.start:matchobj.end]

def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("query", help="Query sequence to search and trim")
    p.add_argument(
        "--input_file", type=argparse.FileType("r"), default=sys.stdin)
    p.add_argument(
        "--trimmed_output_file", type=argparse.FileType("w"), default=sys.stdout)
    p.add_argument(
        "--stats_output_file", type=argparse.FileType("w"), default=sys.stderr)

    # Matching parameters
    p.add_argument(
        "--max_mismatch", type=int, default=0,
        help="Maximum number of mismatches in complete match")
    p.add_argument(
        "--min_partial", type=int, default=0,
        help=(
            "Minimum length of partial sequence match. "
            "Skip partial matching if 0."))
    p.add_argument(
        "--min_pct_id", type=float, default=70.0,
        help="Minimum percent identity in alignment stage")

    # Overall program behavior
    p.add_argument(
        "--min_trimmed_length", type=int, default=1,
        help="Minimum length of trimmed output sequences")
    p.add_argument(
        "--alignment_stages", type=int, default=1,
        help= "Number of pairwise alignment stages")
    p.add_argument(
        "--alignment_dir",
        help=(
            "Directory for files in alignment stage.  If not provided, "
            "a temporary directory will be used."))
    p.add_argument(
        "--cores", type=int, default = 1,
        help="Number of CPU cores to use in alignment stage")
    p.add_argument(
        "--reverse_complement_query", action="store_true",
        help="Reverse complement the query seq before search")
    p.add_argument(
        "--trim_right", action="store_true",
        help="Trim right side, rather than left side of sequences")

    args = p.parse_args(argv)

    seqs = TrimmableSeqs.from_fasta(args.input_file)
    writer = Writer(args.trimmed_output_file, args.stats_output_file)
    app = TrimraggedApp(seqs, args.trim_right, writer, args.min_trimmed_length)

    queryset = deambiguate(args.query)
    if args.reverse_complement_query:
        queryset = [reverse_complement(q) for q in queryset]

    if args.alignment_dir is None:
        alignment_dir = tempfile.mkdtemp()
    else:
        alignment_dir = args.alignment_dir
        if os.path.exists(alignment_dir):
            if not os.path.isdir(alignment_dir):
                raise RuntimeError(
                    "{0} exists and is not a directory".format(alignment_dir))
        else:
            os.mkdir(alignment_dir)

    app.matchers.append(CompleteMatcher(queryset, args.max_mismatch))
    # TODO: Partial matching does not work on right hand side
    if args.min_partial > 0:
        app.matchers.append(PartialMatcher(queryset, args.min_partial))
    for n in range(args.alignment_stages):
        app.matchers.append(AlignmentMatcher(
            alignment_dir,
            min_pct_id = args.min_pct_id,
            cores = args.cores,
            suffix = str(n),
        ))

    app.run()

    if args.alignment_dir is None:
        shutil.rmtree(alignment_dir)

AMBIGUOUS_BASES = {
    "T": "T",
    "C": "C",
    "A": "A",
    "G": "G",
    "R": "AG",
    "Y": "TC",
    "M": "CA",
    "K": "TG",
    "S": "CG",
    "W": "TA",
    "H": "TCA",
    "B": "TCG",
    "V": "CAG",
    "D": "TAG",
    "N": "TCAG",
    }


# Ambiguous base codes for all bases EXCEPT the key
AMBIGUOUS_BASES_COMPLEMENT = {
    "T": "V",
    "C": "D",
    "A": "B",
    "G": "H",
    }


def deambiguate(seq):
    nt_choices = [AMBIGUOUS_BASES[x] for x in seq]
    return ["".join(c) for c in itertools.product(*nt_choices)]


COMPLEMENT_BASES = {
    "T": "A",
    "C": "G",
    "A": "T",
    "G": "C",
    }


def reverse_complement(seq):
    rc = [COMPLEMENT_BASES[x] for x in seq]
    rc.reverse()
    return ''.join(rc)
