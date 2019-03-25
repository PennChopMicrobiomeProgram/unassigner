from unassign.util import count_while_equal, count_matching_pairs


class AlignedSubjectQuery(object):
    def __init__(self, qseq, sseq):
        self.query_id, self.query_seq, self.query_len = qseq
        self.subject_id, self.subject_seq, self.subject_len = sseq
        assert(len(self.query_seq) == len(self.subject_seq))

    def pairs_query(self, start_idx = 0, end_idx = None):
        if end_idx is None:
            end_idx = self.query_len
        query_idx = 0
        for q, s in zip(self.query_seq, self.subject_seq):
            if query_idx >= end_idx:
                break
            if query_idx >= start_idx:
                yield (q, s)
            if q != "-":
                query_idx += 1

    def count_matches(self):
        """Count regional and total matches in an alignment. 
        Parameters 
        ----------
        start : start position in query sequence
        end : end position in query sequence
        Returns
        -------
        tuple containing two ints:
        Number of matching positions and 
        total number of query nucleotides in the alignment.        
        """
        matches, total_query, _ = count_matching_pairs(
            zip(self.query_seq, self.subject_seq))
        return matches, self.query_len
