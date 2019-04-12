import itertools
import operator

class AlignedPair(object):
    def __init__(self, qseq, sseq):
        self.query_id, self.query_seq = qseq
        self.subject_id, self.subject_seq = sseq
        assert(len(self.query_seq) == len(self.subject_seq))

    @property
    def alignment_len(self):
        return len(self.query_seq)

    @property
    def query_len(self):
        return len(self.query_seq) - self.query_seq.count("-")

    @property
    def subject_len(self):
        return len(self.subject_seq) - self.subject_seq.count("-")

    @property
    def unaligned_query_seq(self):
        return self.query_seq.replace("-", "")

    @property
    def unaligned_subject_seq(self):
        return self.subject_seq.replace("-", "")

    def count_matches(self):
        return sum(q == s for q, s in zip(self.query_seq, self.subject_seq))

class AlignedRegion:
    def __init__(self, alignment, start_idx, end_idx):
        assert(start_idx >= 0)
        assert(start_idx <= alignment.alignment_len)
        assert(end_idx >= 0)
        assert(end_idx <= alignment.alignment_len)
        assert(start_idx <= end_idx)
        self.alignment = alignment
        self.start_idx = start_idx
        self.end_idx = end_idx

    def trim_ends(self):
        qseq = self.alignment.query_seq[self.start_idx:self.end_idx]
        sseq = self.alignment.subject_seq[self.start_idx:self.end_idx]
        return AlignedPair(
            (self.alignment.query_id, qseq),
            (self.alignment.subject_id, sseq))

    def trim_left(self, include_region=False):
        if include_region:
            idx = self.start_idx
        else:
            idx = self.end_idx
        qseq = self.alignment.query_seq[idx:]
        sseq = self.alignment.subject_seq[idx:]
        return AlignedPair(
            (self.alignment.query_id, qseq),
            (self.alignment.subject_id, sseq))

    def trim_right(self, include_region=False):
        if include_region:
            idx = self.end_idx
        else:
            idx = self.start_idx
        qseq = self.alignment.query_seq[:idx]
        sseq = self.alignment.subject_seq[:idx]
        return AlignedPair(
            (self.alignment.query_id, qseq),
            (self.alignment.subject_id, sseq))

    def in_alignment(self):
        return (self.start_idx, self.end_idx)

    def in_subject(self):
        return in_aligned_seq(
            self.alignment.subject_seq, self.start_idx, self.end_idx)

    def in_query(self):
        return in_aligned_seq(
            self.alignment.query_seq, self.start_idx, self.end_idx)

    def subject_offset(self):
        qstart, qend = self.in_subject()
        qlen = self.alignment.subject_len
        if (qstart == 0) and (qend == 0):
            # Subject is off to right side, offset is positive
            subject_right = self.alignment.subject_seq[self.end_idx:]
            return count_endgaps(subject_right)
        if (qstart == qlen) and (qend == qlen):
            # Subject is off to left side, offset is negative
            subject_left = self.alignment.subject_seq[:self.start_idx]
            return -count_endgaps(subject_left, reverse=True)
        return 0

    def query_offset(self):
        qstart, qend = self.in_query()
        qlen = self.alignment.query_len
        if (qstart == 0) and (qend == 0):
            # Query is off to right side, offset is positive
            query_right = self.alignment.query_seq[self.end_idx:]
            return count_endgaps(query_right)
        if (qstart == qlen) and (qend == qlen):
            # Query is off to left side, offset is negative
            query_left = self.alignment.query_seq[:self.start_idx]
            return -count_endgaps(query_left, reverse=True)
        return 0

    @classmethod
    def without_endgaps(cls, a):
        left_endgaps = max(
            count_endgaps(a.query_seq),
            count_endgaps(a.subject_seq))
        right_endgaps = max(
            count_endgaps(a.query_seq, reverse=True),
            count_endgaps(a.subject_seq, reverse=True))
        start_idx = left_endgaps
        end_idx = a.alignment_len - right_endgaps
        return cls(a, start_idx, end_idx)

    @classmethod
    def from_subject(cls, a, subject_start_idx=0, subject_end_idx=None):
        if subject_end_idx is None:
            subject_end_idx = a.subject_len
        assert(subject_start_idx >= 0)
        assert(subject_start_idx <= a.subject_len)
        assert(subject_end_idx >= 0)
        assert(subject_end_idx <= a.subject_len)
        assert(subject_start_idx <= subject_end_idx)
        start_idx = aligned_start_idx(a.subject_seq, subject_start_idx)
        end_idx = aligned_end_idx(a.subject_seq, subject_end_idx)
        return cls(a, start_idx, end_idx)

    @classmethod
    def from_query(cls, a, query_start_idx=0, query_end_idx=None):
        if query_end_idx is None:
            query_end_idx = a.query_len
        assert(query_start_idx >= 0)
        assert(query_start_idx <= a.query_len)
        assert(query_end_idx >= 0)
        assert(query_end_idx <= a.query_len)
        assert(query_start_idx <= query_end_idx)
        start_idx = aligned_start_idx(a.query_seq, query_start_idx)
        end_idx = aligned_end_idx(a.query_seq, query_end_idx)
        return cls(a, start_idx, end_idx)


def in_aligned_seq(seq, start_idx, end_idx):
    seq_start = seq[:start_idx]
    seq_start_idx = start_idx - seq_start.count("-")
    seq_end = seq[:end_idx]
    seq_end_idx = end_idx - seq_end.count("-")
    return seq_start_idx, seq_end_idx

def aligned_end_idx(seq, end_idx):
    # Work in reverse mode, find the complimentary start_idx
    ungapped_seq_len = len(seq) - seq.count("-")
    rev_start_idx = ungapped_seq_len - end_idx
    rev_seq = seq[::-1]
    rev_alignment_start_idx = aligned_start_idx(rev_seq, rev_start_idx)
    return len(seq) - rev_alignment_start_idx

def aligned_start_idx(seq, start_idx):
    ungapped_seq_len = len(seq) - seq.count("-")
    # Special case: start_idx is equal to ungapped sequence length
    # We will never count enough bases to reach the start_idx
    # Reverse the sequence and find the start position for index 0,
    # then subtract from the alignment length to get the forward-facing
    # index after the last base.
    if start_idx == ungapped_seq_len:
        rev_seq = seq[::-1]
        rev_alignment_start_idx = aligned_start_idx(rev_seq, 0)
        return len(seq) - rev_alignment_start_idx
    start_idxs = [n - 1 for n in cumulative_bases(seq)]
    return start_idxs.index(start_idx)

def cumulative_bases(seq):
    return itertools.accumulate(int(b != "-") for b in seq)

def iter_len(xs):
    n = 0
    for _ in xs:
        n += 1
    return n

def count_endgaps(seq, reverse = False):
    if reverse:
        seq = reversed(seq)
    return iter_len(itertools.takewhile(lambda x: x == "-", seq))

