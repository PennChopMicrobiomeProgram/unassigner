from unassigner.alignment import AlignedPair

class HitExtender:
    def __init__(self, query_seqs, ref_seqs):
        self.query_seqs = dict(query_seqs)
        self.ref_seqs = dict(ref_seqs)

    def extend_hit(self, hit):
        # Handle the simple case where the local alignment covers both
        # sequences completely
        if self._is_global(hit):
            return AlignedPair(
                (hit['qseqid'], hit['qseq']),
                (hit['sseqid'], hit['sseq']))

        # We are going to need some repair or realignment.
        qseq = self.query_seqs[hit['qseqid']]
        assert(len(qseq) == hit['qlen'])
        sseq = self.ref_seqs[hit['sseqid']]
        assert(len(sseq) == hit['slen'])

        # We expect this always to be true for the vsearch aligner.
        # In the past, we had to re-do the alignment ourselves if it
        # didn't come back from the external program as a semi-global
        # alignment. I'm leaving this in place to double-check that we
        # have a semi-global alignment. If we go around adding endgaps
        # to an alignment that's not semi-global, we are sure to
        # generate some difficult-to-diagnose bugs.
        assert self._is_semiglobal(hit)

        qleft, sleft = self._add_endgaps_left(hit, qseq, sseq)
        qright, sright = self._add_endgaps_right(hit, qseq, sseq)
        aligned_qseq = qleft + hit['qseq'] + qright
        aligned_sseq = sleft + hit['sseq'] + sright
        return AlignedPair(
                (hit['qseqid'], aligned_qseq),
                (hit['sseqid'], aligned_sseq))

    @staticmethod
    def _is_global(hit):
        return (
            (hit['qstart'] == 1) and \
            (hit['sstart'] == 1) and \
            (hit['qend'] == hit['qlen']) and \
            (hit['send'] == hit['slen']))

    @staticmethod
    def _is_semiglobal(hit):
        more_to_the_left = (hit['qstart'] > 1) and \
                           (hit['sstart'] > 1)
        more_to_the_right = (hit['qend'] < hit['qlen']) and \
                            (hit['send'] < hit['slen'])
        return not (more_to_the_left or more_to_the_right)

    @staticmethod
    def _add_endgaps_left(hit, qseq, sseq):
        # No repair needed
        if (hit['qstart'] == 1) and (hit['sstart'] == 1):
            return ("", "")
        # Query hanging off to the left
        if (hit['qstart'] > 1) and (hit['sstart'] == 1):
            endgap_len = hit['qstart'] - 1
            return (qseq[:endgap_len], "-" * endgap_len)
        # Subject hanging off to the left
        if (hit['qstart'] == 1) and (hit['sstart'] > 1):
            endgap_len = hit['sstart'] - 1
            return ("-" * endgap_len, sseq[:endgap_len])
        # Anything not meeting these conditions is bad
        if (hit['qstart'] > 1) and (hit['sstart'] > 1):
            raise ValueError("Unaligned sequence on left")
        raise ValueError("Query or subject start position less than 1")

    @staticmethod
    def _add_endgaps_right(hit, qseq, sseq):
        # No repair needed
        if (hit['qend'] == hit['qlen']) and (hit['send'] == hit['slen']):
            return ("", "")
        # Query hanging off to the right
        if (hit['qend'] < hit['qlen']) and (hit['send'] == hit['slen']):
            endgap_len = hit['qlen'] - hit['qend']
            return (qseq[-endgap_len:], "-" * endgap_len)
        # Subject hanging off to the right
        if (hit['qend'] == hit['qlen']) and (hit['send'] < hit['slen']):
            endgap_len = hit['slen'] - hit['send']
            return ("-" * endgap_len, sseq[-endgap_len:])
        # Anything not meeting these conditions is bad
        if (hit['qend'] < hit['qlen']) and (hit['send'] < hit['qlen']):
            raise ValueError("Unaligned sequence on right")
        raise ValueError("Query or subject end position greater than length")
