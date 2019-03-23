from unassign.util import count_while_equal, count_matching_pairs


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

    @classmethod
    def from_blast_hit(cls, hit):
        return cls(
            (hit['qseqid'], hit['qseq'], hit['qlen']),
            (hit['sseqid'], hit['sseq'], hit['slen']))

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
        
        matches, total_query, _ = count_matching_pairs(zip(self.query_seq, self.subject_seq))
        return matches, len(self.query_seq)
