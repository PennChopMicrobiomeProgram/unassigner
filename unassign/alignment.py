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

    def region_subject_to_query(self, subject_start_idx=0, subject_end_idx=None):
        "The region of the query inside the specified region of the subject."
        if subject_end_idx is None:
            subject_end_idx = self.subject_len
        query_idx = 0
        subject_idx = 0
        query_start_idx = None
        for q, s in zip(self.query_seq, self.subject_seq):
            if (query_start_idx is None) and \
               (subject_idx >= subject_start_idx) and \
               (s != "-"):
                query_start_idx = query_idx
            if s != "-":
                subject_idx += 1
            if q != "-":
                query_idx += 1
            if subject_idx >= subject_end_idx:
                break
        return query_start_idx, query_idx

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

    def region_query_to_subject(self, query_start_idx=0, query_end_idx=None):
        "The region of the query inside the specified region of the query."
        if query_end_idx is None:
            query_end_idx = self.query_len
        query_idx = 0
        subject_idx = 0
        subject_start_idx = None
        for q, s in zip(self.query_seq, self.subject_seq):
            if (subject_start_idx is None) and \
               (query_idx >= query_start_idx) and \
               (q != "-"):
                subject_start_idx = subject_idx
            if s != "-":
                subject_idx += 1
            if q != "-":
                query_idx += 1
            if query_idx >= query_end_idx:
                break
        return subject_start_idx, subject_idx

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
