from unassign.util import count_while_equal


class Alignment(object):
    def __init__(self, qseq, sseq):
        self.query_id, self.query_seq, self.query_len = qseq
        self.subject_id, self.subject_seq, self.subject_len = sseq
        assert(len(self.query_seq) == len(self.subject_seq))

    @staticmethod
    def start_idx(subject_seq, query_seq):
        return max(count_while_equal(query_seq, "-"),
                   count_while_equal(subject_seq, "-"))

    @staticmethod
    def end_idx(subject_seq, query_seq):
        return len(query_seq) - max(
            count_while_equal(reversed(query_seq), "-"),
            count_while_equal(reversed(subject_seq), "-"))

    def get_pairs(self, start=None, end=None):
        astart = self.start_idx
        if (start is None) or (start < astart):
            start = astart

        aend = self.end_idx
        if (end is None) or (end > aend):
            end = aend

        return list(zip(self.query_seq, self.subject_seq))[start:end]

    def count_matches(self, start=None, end=None):
        """Count regional and total matches in BLAST hit. 
        Parameters 
        ----------
        start : start position in query sequence
        end : end position in query sequence

        Returns
        -------
        tuple containing two ints:
        Number of matching positions and 
        total number of query nucleotides in the alignment.
        
        Notes
        -----
        Because sequence positions are indexed from 1 in BLAST, the
        start and end positions are indexed starting from 1.
        """
        
        if start is None:
            start = 1
        if end is None:
            end = len(self.query_seq)

        total = 0
        matches = 0

        qpos = start
        for qchar, hchar in zip(self.query_seq, self.subject_seq):
            if start <= qpos <= end: ## TODO: This can be simplified if we are sure we are given the full alignment for the query
                total += 1
                if qchar == hchar:
                    matches += 1
            if qchar != '-':
                qpos += 1

        return matches, total
