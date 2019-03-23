import abc
import argparse
import collections
import itertools
import os
import subprocess
import sys

from unassign.parse import parse_fasta, write_fasta

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

    def get_matched_recs(self):
        for rep_seq_id in sorted(self.seq_ids):
            if rep_seq_id in self.matches:
                seq = self.seqs[rep_seq_id]
                yield rep_seq_id, seq

    def get_unmatched_recs(self):
        for rep_seq_id in sorted(self.seq_ids):
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
    "PrimerMatch", ["start", "end", "message"])


class Matcher(abc.ABC):
    def __init__(self, queryset):
        self.queryset = queryset

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
                    return PrimerMatch(start_idx, end_idx, msg)


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
                return PrimerMatch(0, end_idx, "Partial")

class AlignmentMatcher(Matcher):
    def __init__(self, min_pct_id = 90, keep_files = False, cores = 1):
        self.min_pct_id = min_pct_id
        self.keep_files = keep_files
        self.cores = cores

    def find_in_seqs(self, seqs):
        if seqs.all_matched():
            raise StopIteration()

        # Write query
        query_fp = ".trimragged.query.fa"
        with open(query_fp, "w") as f:
            write_fasta(f, seqs.get_unmatched_recs())

        # Write database
        database_fp = ".trimragged.database.fa"
        with open(database_fp, "w") as f:
            write_fasta(f, seqs.get_matched_recs())
        subprocess.check_call(
            ["makeblastdb", "-dbtype", "nucl", "-in", database_fp],
            stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL,
        )

        results_fp = ".trimragged.results.txt"
        command = [
            "blastn",
            "-outfmt", "7",
            "-query", query_fp,
            "-db", database_fp,
            "-out", results_fp,
            "-max_target_seqs", "1",
        ]
        if self.cores > 1:
            command.extend(["-num_threads", str(self.cores)])
        subprocess.check_call(command)

        # Read output file, determine matches
        with open(results_fp, "r") as f:
            for hit in parse_blast7(f):
                if hit["pident"] < self.min_pct_id:
                    continue
                subj_id = hit["sseqid"]
                subj_match = seqs.matches[subj_id]

                query_alignment_start_idx = hit["qstart"] - 1
                subj_alignment_start_idx = hit["sstart"] - 1

                query_start_idx = subj_match.start - \
                    (subj_alignment_start_idx - query_alignment_start_idx)
                if query_start_idx < 0:
                    query_start_idx = 0

                query_end_idx = subj_match.end - \
                    (subj_alignment_start_idx - query_alignment_start_idx)
                matchobj = PrimerMatch(
                    query_start_idx, query_end_idx, "Alignment")
                query_id = hit["qseqid"]
                yield query_id, matchobj

        if not self.keep_files:
            os.remove(database_fp)
            os.remove(query_fp)
            os.remove(results_fp)


def parse_blast7(f):
    BLAST_FIELDS = [
        "qseqid", "sseqid", "pident", "length",
        "mismatch", "gapopen", "qstart", "qend",
        "sstart", "send", "evalue", "bitscore",
    ]
    BLAST_FIELD_TYPES = [
        str, str, float, int,
        int, int, int, int,
        int, int, float, float,
    ]
    for line in f:
        line = line.strip()
        if line.startswith("#"):
            continue
        vals = line.split("\t")
        vals = [fn(v) for fn, v in zip(BLAST_FIELD_TYPES, vals)]
        yield dict(zip(BLAST_FIELDS, vals))


def is_digit(char):
    return char in "1234567890"


def pairs(xs):
    args = [iter(xs)] * 2
    return zip(*args)


def partial_seqs(seq, min_length):
    max_start_idx = len(seq) - min_length + 1
    for start_idx in range(1, max_start_idx):
        yield seq[start_idx:]


class TrimraggedApp(object):
    def __init__(self, seqs, trim_fcn, writer):
        self.seqs = seqs
        self.trim_fcn = trim_fcn
        self.writer = writer

    def apply_matcher(self, matcher):
        matches = matcher.find_in_seqs(self.seqs)
        for rep_seq_id, matchobj in matches:
            if matchobj is not None:
                self.seqs.register_match(rep_seq_id, matchobj)
                for seq_id, seq in self.seqs.get_replicate_recs(rep_seq_id):
                    desc = self.seqs.get_desc(seq_id)
                    trimmed_seq = self.trim_fcn(seq, matchobj)
                    self.writer.write_trimmed(desc, trimmed_seq)
                    self.writer.write_stats(seq_id, matchobj)

    def finish(self):
        for rep_seq_id, seq in self.seqs.get_unmatched_recs():
            for seq_id in self.seqs.get_replicate_ids(rep_seq_id):
                desc = self.seqs.get_desc(seq_id)
                self.writer.write_untrimmed(desc, seq)
                self.writer.write_stats(seq_id, None)


class Writer(object):
    def __init__(self, trimmed_file, stats_file, untrimmed_file):
        self.trimmed_file = trimmed_file
        self.stats_file = stats_file
        self.untrimmed_file = untrimmed_file

    def write_trimmed(self, desc, seq):
        self.trimmed_file.write(">{0}\n{1}\n".format(desc, seq))

    def write_stats(self, seq_id, matchobj):
        if matchobj is None:
            self.stats_file.write("{0}\tUnmatched\tNA\tNA\n".format(seq_id))
        else:
            self.stats_file.write("{0}\t{1}\t{2}\t{3}\n".format(
                seq_id, matchobj.message, matchobj.start, matchobj.end))

    def write_untrimmed(self, desc, seq):
        if self.untrimmed_file is not None:
            self.untrimmed_file.write(">{0}\n{1}\n".format(desc, seq))


def trim_left(seq, matchobj):
    return seq[matchobj.end:]

def trim_right(seq, matchobj):
    return seq[:matchobj.start]

def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument(
        "--input_file", type=argparse.FileType("r"), default=sys.stdin)
    p.add_argument(
        "--trimmed_output_file", type=argparse.FileType("w"), default=sys.stdout)
    p.add_argument(
        "--stats_output_file", type=argparse.FileType("w"), default=sys.stderr)
    p.add_argument(
        "--unmatched_output_file", type=argparse.FileType("w"))
    p.add_argument("--query", required=True)

    # Parameters for each step
    p.add_argument(
        "--max_mismatch", type=int, default=2,
        help="Maximum number of mismatches in complete match")
    p.add_argument(
        "--min_partial", type=int, default=5,
        help="Minimum length of partial sequence match")
    p.add_argument(
        "--min_pct_id", type=float, default=80.0,
        help="Minimum percent identity in alignment stage")

    # Overall program behavior
    p.add_argument(
        "--skip_alignment", action="store_true",
        help= "Skip pairwise alignment stage")
    p.add_argument(
        "--skip_partial", action="store_true",
        help="Skip partial alignment stage")
    p.add_argument(
        "--keep_alignment_files", action="store_true",
        help="Keep database, query, and results files from alignment stage")
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
    writer = Writer(
        args.trimmed_output_file, args.stats_output_file,
        args.unmatched_output_file)
    if args.trim_right:
        trim_fcn = trim_right
    else:
        trim_fcn = trim_left
    app = TrimraggedApp(seqs, trim_fcn, writer)

    queryset = deambiguate(args.query)
    if args.reverse_complement_query:
        queryset = [reverse_complement(q) for q in queryset]

    matchers = [CompleteMatcher(queryset, args.max_mismatch)]
    if not args.skip_partial:
        matchers.append(PartialMatcher(queryset, args.min_partial))
    if not args.skip_alignment:
        matchers.append(
            AlignmentMatcher(
                args.min_pct_id, args.keep_alignment_files,
                args.cores))

    for m in matchers:
        app.apply_matcher(m)
    app.finish()


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
