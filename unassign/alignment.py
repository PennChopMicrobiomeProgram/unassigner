class AlignedSubjectQuery(object):
    def __init__(self, qseq, sseq):
        self.query_id, self.query_seq = qseq
        self.subject_id, self.subject_seq = sseq
        assert(len(self.query_seq) == len(self.subject_seq))

    @property
    def query_len(self):
        return len(self.query_seq) - self.query_seq.count("-")

    @property
    def subject_len(self):
        return len(self.subject_seq) - self.subject_seq.count("-")

    def pairs_query(self, start_idx = 0, end_idx = None):
        if end_idx is None:
            end_idx = self.query_len
        query_idx = 0
        for q, s in zip(self.query_seq, self.subject_seq):
            if q != "-":
                query_idx += 1
            if (query_idx == end_idx) and (q == "-"):
                break
            if query_idx > end_idx:
                break
            if query_idx > start_idx:
                yield (q, s)

    def pairs_subject(self, start_idx = 0, end_idx = None):
        if end_idx is None:
            end_idx = self.subject_len
        subject_idx = 0
        for q, s in zip(self.query_seq, self.subject_seq):
            if s != "-":
                subject_idx += 1
            if (subject_idx == end_idx) and (s == "-"):
                break
            if subject_idx > end_idx:
                break
            if subject_idx > start_idx:
                yield (q, s)

def exact_match(a, b):
    return a == b

def count_matches(sequence_pairs, match_fcn=exact_match):
    match = 0
    total_a = 0
    total_b = 0
    for a, b in sequence_pairs:
        if match_fcn(a, b):
            match += 1
        if a != "-":
            total_a += 1
        if b != "-":
            total_b += 1
    return match, total_a, total_b
