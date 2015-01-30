from unassign.util import count_while_equal


class Alignment(object):
    def __init__(self, qseq, sseq):
        self.query_id, self.query_seq = qseq
        self.subject_id, self.subject_seq = sseq
        assert(len(self.query_seq) == len(self.subject_seq))

    @property
    def start_idx(self):
        return max(
            count_while_equal(self.query_seq, "-"),
            count_while_equal(self.subject_seq, "-"))

    @property
    def end_idx(self):
        return len(self.query_seq) - max(
            count_while_equal(reversed(self.query_seq), "-"),
            count_while_equal(reversed(self.subject_seq), "-"))

    def get_pairs(self, start=None, end=None):
        astart = self.start_idx
        if (start is None) or (start < astart):
            start = astart

        aend = self.end_idx
        if (end is None) or (end > aend):
            end = aend

        return zip(self.query_seq, self.subject_seq)[start:end]
