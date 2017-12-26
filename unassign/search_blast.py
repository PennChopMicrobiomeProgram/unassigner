from __future__ import division
import itertools
import subprocess
import tempfile
from Bio import pairwise2

from unassign.parse import write_fasta, load_fasta, parse_fasta
from unassign.util import uniq
from unassign.alignment import Alignment

BLAST_FMT = (
    "qseqid sseqid pident length mismatch gapopen "
    "qstart qend sstart send qlen slen qseq sseq")
BLAST_FIELDS = BLAST_FMT.split()
BLAST_FIELD_TYPES = [
    str, str, float, int, int, int,
    int, int, int, int, int, int, str, str]


class BlastAlignment(Alignment):
    def __init__(self, hit):
        self._hit = hit
        super(BlastAlignment, self).__init__(
            (hit['qseqid'], hit['qseq'], hit['qlen']), (hit['sseqid'], hit['sseq'], hit['slen']))

    def get_local_pairs(self):
        return zip(self._hit['qseq'], self._hit['sseq'])
    
class SemiGlobalAlignment(Alignment):
    def __init__(self, query_id, qseq_orj, subject_id, sseq_orj):
        query_seq, subject_seq, qlen, slen  = self._get_aligned_seqs(qseq_orj, sseq_orj)
        super(SemiGlobalAlignment, self).__init__(
            (query_id, query_seq, qlen), (subject_id, subject_seq, slen))

    def _get_aligned_seqs(self, qseq_orj, sseq_orj):
        alignment = pairwise2.align.globalms(sseq_orj, qseq_orj,
                                             5, -4, -10, -0.5, #match, mismatch, gapopen, gapextend #### TODO: make these configurable
                                             penalize_end_gaps=False, one_alignment_only=True)
        subj_seq, query_seq = self._trim_global_alignment(alignment[0][0], alignment[0][1])
        return query_seq, subj_seq, len(qseq_orj), len(sseq_orj)

    def _trim_global_alignment(self, subj_seq, query_seq):
        # trim only the gaps on either side of the QUERY sequence. If the subject has gaps
        # at the beginning or the end they will be counted as mismatches!!
        non_indel_indices = [i for i, c in enumerate(query_seq) if c!='-']
        triml = non_indel_indices[0]
        trimr = min(non_indel_indices[-1] + 1, len(query_seq))
        return subj_seq[triml:trimr], query_seq[triml:trimr]

class BlastAligner(object):
    """Align sequences with BLAST."""

    def __init__(self, species_fp):
        self.species_fp = species_fp
        self.species_max_hits = 1
        self.species_input_fp = None
        self.species_output_fp = None

        self.num_cpus = 1 #### TODO: make this configurable

    def search_species(self, seqs):
        """Search species typestrains for match to query sequences."""
        return self._search(
            seqs, self.species_fp, self.species_max_hits,
            self.species_input_fp, self.species_output_fp)

    def _get_species_seqs(self, hits):
        """Fetch seqs for each species in the list of hits.

        Each unique species in the list is returned only once in the
        list of results, but the order of hits is preserved to
        facilitate caching.
        """
        species_ids = uniq(x.subject_id for x in hits)
        seqs = load_fasta(self.species_fp)
        return [(x, seqs[x]) for x in species_ids]

    def _search(self, seqs, db, max_hits, input_fp, output_fp):
        if input_fp is None:
            infile = tempfile.NamedTemporaryFile(mode="w+t", encoding="utf-8")
            write_fasta(infile, seqs)
            infile.seek(0)
            input_fp = infile.name
        else:
            with open(input_fp, "w") as f:
                write_fasta(f, seqs)

        if output_fp is None:
            outfile = tempfile.NamedTemporaryFile()
            output_fp = outfile.name

        self._call(
            input_fp, db, output_fp,
            max_target_seqs=max_hits,
            num_threads=self.num_cpus)
        return self._load(output_fp, seqs, db)

    @classmethod
    def _load(self, output_fp, seqs, db):
        """Load hits from an output file."""
        with open(output_fp) as f:                
            hits = [self._polish_alignment(hit, seqs, db) for hit in self._parse(f)]
        return hits

    @classmethod
    def _parse(self, f):
        """Parse a BLAST output file."""
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            vals = line.split("\t")
            vals = [fn(v) for fn, v in zip(BLAST_FIELD_TYPES, vals)]
            yield dict(zip(BLAST_FIELDS, vals))

    @staticmethod
    def _index(fasta_fp):
        return subprocess.check_call([
            "makeblastdb",
            "-dbtype", "nucl",
            "-in", fasta_fp,
            ])

    @classmethod
    def _polish_alignment(self, hit, seqs, db):
        bool_gaps_left = hit['qstart'] > 1
        bool_gaps_right = hit['qend'] < hit['qlen']
        if bool_gaps_left or bool_gaps_right:
            qseq_orj = self._get_orj_query_seq(seqs, hit['qseqid'])
            sseq_orj = self._get_orj_subject_seq(db, hit['sseqid'])
            return SemiGlobalAlignment(hit['qseqid'], qseq_orj,  hit['sseqid'], sseq_orj)
        else:
            return BlastAlignment(hit)
            
    @classmethod
    def _get_orj_query_seq(self, seqs, query_id):
        for seq in seqs: ##### TODO: change seqs to a dictionary to save time here in read_fasta
            if seq[0] == query_id:
                return seq[1]

    @classmethod
    def _get_orj_subject_seq(self, db, subject_id):
        subject_outfile = tempfile.NamedTemporaryFile()
        subject_outfile_fp = subject_outfile.name
        args = ["blastdbcmd",
                "-db", db,
                "-entry", subject_id,
                "-out", subject_outfile_fp
        ]
        subprocess.check_call(args)
        with open(subject_outfile_fp) as f:
            return list(parse_fasta(f, trim_desc=True))[0][1]
            
    def _call(self, query_fp, database_fp, output_fp, **kwargs):
        """Call the BLAST program."""
        args = [
            "blastn",
            "-evalue", "1e-5",
            "-outfmt", "6 " + BLAST_FMT,
            ]
        for arg, val in kwargs.items():
            arg = "-" + arg
            if val is None:
                args.append(arg)
            else:
                args += [arg, str(val)]
        args += [
            "-query", query_fp,
            "-db", database_fp,
            "-out", output_fp,
            ]
        subprocess.check_call(args)
