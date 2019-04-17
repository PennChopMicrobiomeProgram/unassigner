def species_probability(self, species_alignment, refseq_alignments):
    """Compute probability of species attribution.
    """
    query_id = species_alignment.query_id
    species_id = species_alignment.subject_id
    start = species_alignment.start_pos
    end = species_alignment.end_pos
    for r in refseq_alignments:
        refseq_id = r.subject_id
        a, b = r.count_matches(start, end)
        c, d = r.count_matches()
        yield query_id, species_id, a, b, r.subject_id, c, d
