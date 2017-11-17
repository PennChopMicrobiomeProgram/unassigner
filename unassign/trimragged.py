import argparse
import collections
import itertools
import subprocess
import sys

from unassign.parse import parse_fasta, write_fasta

class SeqRecord(object):
    def __init__(self, desc, seq):
        self.desc = desc
        self.seq = seq

    @property
    def seq_id(self):
        return self.desc.split()[0]

    def write_fasta(self, f):
        f.write(">{0}\n{1}\n".format(self.desc, self.seq))


class QueryMatch(object):
    def __init__(self, rec, match_start, match_end, message):
        self.rec = rec
        self.start = match_start
        self.end = match_end
        self.message = message

    def write_stats(self, f):
        fields = (self.rec.seq_id, self.message, self.start, self.end)
        f.write("\t".join(map(str, fields)))
        f.write("\n")

    def trim_left(self):
        return SeqRecord(self.rec.desc, self.rec.seq[self.end:])


class Matcher(object):
    def __init__(self, queryset):
        self.queryset = queryset

    def partition_matching_seqs(self, recs):
        matched = []
        unmatched = []
        for rec in recs:
            match = self.find_match(rec)
            if match is None:
                unmatched.append(rec)
            else:
                matched.append(match)
        return matched, unmatched
                
    def find_match(self, rec):
        raise NotImplemented()


class CompleteMatcher(Matcher):
    def __init__(self, queryset, max_mismatch=2):
        self.queryset = queryset
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
        assert(n_mismatch in [0, 1, 2])
        for query in self.queryset:
            idx_sets = itertools.combinations(range(len(query)), n_mismatch)
            for idx_set in idx_sets:
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

    def find_match(self, rec):
        for n_mismatches, queryset in enumerate(self.mismatched_queryset):
            for query in queryset:
                start_idx = rec.seq.find(query)
                if start_idx > -1:
                    if n_mismatches == 0:
                        msg = "Exact"
                    elif n_mismatches == 1:
                        msg = "Complete, 1 mismatch"
                    else:
                        msg = "Complete, {0} mismatches".format(n_mismatches)
                    return QueryMatch(
                        rec, start_idx, start_idx + len(query), msg)


class PartialMatcher(Matcher):
    def __init__(self, queryset, min_length):
        super().__init__(queryset)
        self.min_length = min_length
        
        self.partial_queries = set()
        for query in self.queryset:
            for partial_query in partial_seqs(query, self.min_length):
                self.partial_queries.add(partial_query)
        
    def find_match(self, rec):
        for partial_query in self.partial_queries:
            if rec.seq.startswith(partial_query):
                return QueryMatch(
                    rec, 0, len(partial_query), "Partial")

class AlignmentMatcher(Matcher):
    def __init__(self, queryset, prev_matches, min_pct_id=80.0):
        super().__init__(queryset)
        self.prev_matches = dict((m.rec.seq_id, m) for m in prev_matches)
        self.min_pct_id = min_pct_id

    def partition_matching_seqs(self, recs):
        # Store recs in dict
        recs = collections.OrderedDict((r.seq_id, r) for r in recs)

        # Write database
        database_fp = "ggdatabase.fa"
        with open(database_fp, "w") as f:
            for m in self.prev_matches.values():
                m.rec.write_fasta(f)

        # Write query
        query_fp = "ggquery.fa"
        with open(query_fp, "w") as f:
            for rec in recs.values():
                rec.write_fasta(f)

        results_fp = "ggresults.txt"
                
        # Run GGsearch
        command = [
            "ggsearch36", "-b", "1", "-d", "1", "-m", "8CB", "-n", 
            query_fp, database_fp]
        with open(results_fp, "w") as f:
            subprocess.check_call(command, stdout=f)
        
        # Read output file, determine matches
        matched = []
        with open(results_fp, "r") as f:
            for query_id, subject_id, pct_id, btop in parse_ggsearch_8CB(f):
                if pct_id < self.min_pct_id:
                    continue
                subject_match = self.prev_matches[subject_id]
                query_start_idx = get_query_idx(btop, subject_match.start)
                query_end_idx = get_query_idx(btop, subject_match.end)
                query_rec = recs[query_id]
                match = QueryMatch(
                    query_rec, query_start_idx, query_end_idx, "Alignment")
                matched.append(match)
                del recs[query_id]

        return matched, list(recs.values())


def get_query_idx(btop, subject_idx):
    for btop_query_idx, btop_subject_idx in btop_idx(btop):
        if btop_subject_idx == subject_idx:
            return btop_query_idx

def is_digit(char):
    return char in "1234567890"

def btop_idx(btop):
    query_idx = 0
    subject_idx = 0
    for is_digit_group, vals in itertools.groupby(btop, is_digit):
        val = "".join(vals)
        if is_digit_group:
            num_matches = int(val)
            for _ in range(num_matches):
                query_idx += 1
                subject_idx += 1
                yield query_idx, subject_idx
        else:
            for qchar, schar in pairs(val):
                if qchar != "-":
                    query_idx += 1
                if schar != "-":
                    subject_idx += 1
                yield query_idx, subject_idx

def pairs(xs):
    xs = iter(xs)
    while xs:
        yield next(xs), next(xs)

def parse_ggsearch_8CB(f):
    for line in f:
        if line.startswith("#"):
            continue
        line = line.rstrip()
        vals = line.split("\t")
        query_id = vals[0]
        subject_id = vals[1]
        pct_id = float(vals[2])
        btop = vals[12]
        # query, subject, btop
        yield query_id, subject_id, pct_id, btop
    
def partial_seqs(seq, min_length):
    max_start_idx = len(seq) - min_length + 1
    for start_idx in range(1, max_start_idx):
        yield seq[start_idx:]


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument(
        "--input_file", type=argparse.FileType("r"), default=sys.stdin)
    p.add_argument(
        "--trimmed_output_file", type=argparse.FileType("w"), default=sys.stdout)
    p.add_argument(
        "--stats_output_file", type=argparse.FileType("w"))
    p.add_argument("--query", required=True)

    p.add_argument(
        "--max_mismatch", type=int, default=2,
        help="Maximum number of mismatches in complete match")
    p.add_argument(
        "--min_partial", type=int, default=5,
        help="Minimum length of partial sequence match"
    )
    args = p.parse_args(argv)

    queryset = deambiguate(args.query)
    
    recs = list(
        SeqRecord(desc, seq) for desc, seq in parse_fasta(args.input_file))
    matched = []
    
    m1 = CompleteMatcher(queryset, args.max_mismatch)
    newly_matched, recs = m1.partition_matching_seqs(recs)
    matched.extend(newly_matched)
    
    m2 = PartialMatcher(queryset, args.min_partial)
    newly_matched, recs = m2.partition_matching_seqs(recs)
    matched.extend(newly_matched)

    m3 = AlignmentMatcher(queryset, matched)
    newly_matched, recs = m3.partition_matching_seqs(recs)
    matched.extend(newly_matched)

    for m in matched:
        trimmed_seq = m.trim_left()
        trimmed_seq.write_fasta(args.trimmed_output_file)

        if args.stats_output_file is not None:
            m.write_stats(args.stats_output_file)

    for rec in recs:
        if args.stats_output_file is not None:
            args.stats_output_file.write(
                "{0}\tUnmatched\tNA\tNA\n".format(rec.seq_id))

        


class BarcodeAssigner(object):
    def __init__(self, samples, mismatches=1, revcomp=True):
        self.samples = samples
        if mismatches not in [0, 1, 2]:
            raise ValueError(
                "Only 0, 1, or 2 mismatches allowed (got %s)" % mismatches)
        self.mismatches = mismatches
        self.revcomp = revcomp
        # Sample names assumed to be unique after validating input data
        self.read_counts = dict((s.name, 0) for s in self.samples)
        self._init_hash()

    def _init_hash(self):
        self._barcodes = {}
        for s in self.samples:
            # Barcodes assumed to be present after validating input data
            if self.revcomp:
                bc = reverse_complement(s.barcode)
            else:
                bc = s.barcode

            # Barcodes assumed to be unique after validating input data
            self._barcodes[bc] = s

            for error_bc in self._error_barcodes(bc):
                # Barcodes not guaranteed to be unique after
                # accounting for errors
                if error_bc in self._barcodes:
                    raise ValueError(
                        "Barcode %s for sample %s matches barcode for "
                        "sample %s with %s mismatches" % (
                            error_bc, self._barcodes[error_bc], s,
                            self.mismatches))
                else:
                    self._barcodes[error_bc] = s

    def _error_barcodes(self, barcode):
        # If the number of mismatches is set to 0, there will be no
        # error barcodes. Immediately stop the iteration.
        if self.mismatches == 0:
            raise StopIteration
        # Each item in idx_sets is a set of indices where mismatches
        # should occur.
        idx_sets = itertools.combinations(range(len(barcode)), self.mismatches)
        for idx_set in idx_sets:
            # Change to list because strings are immutable
            bc = list(barcode)
            # Replace the base at each mismatch position with an
            # ambiguous base specifying all possibilities BUT the one
            # we see.
            for idx in idx_set:
                bc[idx] = AMBIGUOUS_BASES_COMPLEMENT[bc[idx]]
            # Expand to all possibilities for mismatching at this
            # particular set of positions
            for error_bc in deambiguate(bc):
                yield error_bc
        
    def assign(self, seq):
        sample = self._barcodes.get(seq)
        if sample is not None:
            self.read_counts[sample.name] += 1
        return sample


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

        

if __name__ == "__main__":
    main()
