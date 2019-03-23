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


class SemiGlobalAlignment(Alignment):
    def __init__(self, query_id, qseq_orj, subject_id, sseq_orj):
        query_seq, subject_seq  = self._get_aligned_seqs(qseq_orj, sseq_orj)
        super(SemiGlobalAlignment, self).__init__(
            (query_id, query_seq, len(qseq_orj)), (subject_id, subject_seq, len(sseq_orj)))

    def _get_aligned_seqs(self, qseq_orj, sseq_orj):
        alignment = pairwise2.align.globalms(sseq_orj, qseq_orj,
                                             5, -4, -10, -0.5, #match, mismatch, gapopen, gapextend #### TODO: make these configurable
                                             penalize_end_gaps=False, one_alignment_only=True)
        subj_seq, query_seq = self._trim_global_alignment(alignment[0][0], alignment[0][1])
        return query_seq, subj_seq

    def _trim_global_alignment(self, subj_seq, query_seq):
        triml = self.start_idx(subj_seq, query_seq)
        trimr = self.end_idx(subj_seq, query_seq)
        return subj_seq[triml:trimr], query_seq[triml:trimr]

class UnassignAligner(object):
    def __init__(self, species_fp):
        self.species_fp = species_fp
        self.species_max_hits = 1
        self.species_input_fp = None
        self.species_output_fp = None

    def search_species(self, seqs):
        """Search species typestrains for match to query sequences."""
        b = BlastSearch(self.species_fp)
        r = BlastRefiner(seqs, self.species_fp)
        hits = b.search(
            seqs, self.species_max_hits,
            self.species_input_fp, self.species_output_fp)
        for hit in hits:
            yield r.refine_hit(hit)

class BlastSearch:
    def __init__(self, ref_seqs_fp):
        self.ref_seqs_fp = ref_seqs_fp

    def search(self, seqs, max_hits, input_fp=None, output_fp=None):
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
            input_fp, self.ref_seqs_fp, output_fp,
            max_target_seqs=max_hits)

        with open(output_fp) as f:
            for hit in self._parse(f):
                yield hit

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

class BlastRefiner:
    def __init__(self, seqs, db):
        self.seqs = seqs
        self.db = db

    def refine_hit(self, hit):
        bool_gaps_left = hit['qstart'] > 1
        bool_gaps_right = hit['qend'] < hit['qlen']
        if bool_gaps_left or bool_gaps_right:
            qseq_orj = self._get_query_seq(hit['qseqid'])
            sseq_orj = self._get_subject_seq(hit['sseqid'])
            return SemiGlobalAlignment(
                hit['qseqid'], qseq_orj,  hit['sseqid'], sseq_orj)
        else:
            return Alignment.from_blast_hit(hit)

    def _get_query_seq(self, query_id):
        ##### TODO: change seqs to a dictionary to save time here in read_fasta
        for seq in self.seqs:
            if seq[0] == query_id:
                return seq[1]

    def _get_subject_seq(self, subject_id):
        subject_outfile = tempfile.NamedTemporaryFile()
        subject_outfile_fp = subject_outfile.name
        args = ["blastdbcmd",
                "-db", self.db,
                "-entry", subject_id,
                "-out", subject_outfile_fp
        ]
        subprocess.check_call(args)
        with open(subject_outfile_fp) as f:
            return list(parse_fasta(f, trim_desc=True))[0][1]
