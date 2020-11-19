import argparse
import os
import subprocess

import argparse
import collections
import itertools
import operator
import sys
import tempfile

from unassigner.parse import parse_fasta, write_fasta
from unassigner.align import VsearchAligner

class MismatchLocationApp:
    def __init__(
            self, species_file, ref_fp, mismatch_file,
            batch_size=10, num_cpus=None):
        self.typestrain_seqs = list(parse_fasta(species_file, trim_desc=True))
        self.reference_fasta_fp = ref_fp
        self.mismatch_db_file = mismatch_file

        self.min_pct_id = 97.0
        self.num_threads = num_cpus
        self.batch_size = batch_size
        self.max_hits = "10000"

    @property
    def reference_udb_fp(self):
        base_fp, _ = os.path.splitext(self.reference_fasta_fp)
        return base_fp + ".udb"

    def make_reference_udb(self):
        if os.path.exists(self.reference_udb_fp):
            return None
        args = [
            "vsearch",
            "--makeudb_usearch", self.reference_fasta_fp,
            "--output", self.reference_udb_fp,
        ]
        return subprocess.check_call(args)

    def search_reference_seqs(self, query_seqs):
        query_file = tempfile.NamedTemporaryFile(suffix=".fasta", mode="wt")
        write_fasta(query_file, query_seqs)
        query_file.seek(0)

        reference_hits_file = tempfile.NamedTemporaryFile(suffix=".txt", mode="wt")

        # 97.0 --> 0.97
        vsearch_min_id = "{:.2f}".format(self.min_pct_id / 100)
        vsearch_args =[
            "vsearch", "--usearch_global", query_file.name,
            "--db", self.reference_udb_fp,
            "--userout", reference_hits_file.name,
            "--iddef", "2", "--id", vsearch_min_id,
            "--maxaccepts", self.max_hits,
            "--userfields",
            "query+target+id2+alnlen+mism+gaps+qilo+qihi+tilo+tihi+qs+ts+qrow+trow",
        ]
        if self.num_threads:
            vsearch_args.extend(["--threads", str(self.num_threads)])

        subprocess.check_call(vsearch_args)
        reference_hits_file.seek(0)
        return reference_hits_file

    def find_mismathes(self, reference_hits_file):
        ref_hits_file_read = open(reference_hits_file.name)
        hits = VsearchAligner._parse(ref_hits_file_read)
        for hit in hits:
            if hit["pident"] > 97.0:
                query_id = hit["qseqid"]
                subject_id = hit["sseqid"]
                mismatch_pos = mismatch_query_pos(hit)
                yield (query_id, subject_id, mismatch_pos)

    def run(self):
        self.make_reference_udb()
        for query_seqs in group_by_n(self.typestrain_seqs, self.batch_size):
            reference_hits_file = self.search_reference_seqs(query_seqs)
            ref_mismatches = self.find_mismathes(reference_hits_file)
            write_mismatches(self.mismatch_db_file, ref_mismatches)

def group_by_n(xs, n):
    args = [iter(xs)] * n
    for xgroup in itertools.zip_longest(*args, fillvalue=None):
        yield [x for x in xgroup if x is not None]

def write_mismatches(f, ref_mismatches):
    for query_id, subject_id, mismatch_pos in ref_mismatches:
        f.write(query_id)
        f.write("\t")
        f.write(subject_id)
        f.write("\t")
        f.write("\t".join(map(str, mismatch_pos)))
        f.write("\n")


AMBIGUOUS_BASES = {
    "-": "-",
    "T": "T",
    "C": "C",
    "A": "A",
    "G": "G",
    "R": "AG",
    "Y": "TC",
    "M": "CA",
    "K": "TG",
    "S": "CG",
    "W": "TA",
    "H": "TCA",
    "B": "TCG",
    "V": "CAG",
    "D": "TAG",
    "N": "TCAG",
}

def nucleotides_compatible(nt1, nt2):
    if (nt1 == "N") or (nt2 == "N"):
        return True
    a1 = AMBIGUOUS_BASES[nt1]
    a2 = AMBIGUOUS_BASES[nt2]
    return bool(set(a1).intersection(a2))

def hit_matches_by_alignment_pos(hit):
    qseq = hit["qseq"].upper()
    sseq = hit["sseq"].upper()
    query_pos = hit["qstart"]
    for q, s in zip(qseq, sseq):
        is_match = True
        if (q != s) and not nucleotides_compatible(q, s):
            is_match = False
        yield (query_pos, is_match)
        if q != "-":
            query_pos += 1

def hit_matches_by_query_pos(hit):
    query_pos_groups = itertools.groupby(
        hit_matches_by_alignment_pos(hit), operator.itemgetter(0))
    for qpos, matches in query_pos_groups:
        all_matched = all(is_match for _, is_match in matches)
        yield qpos, all_matched

def mismatch_query_pos(hit):
    for qpos, is_match in hit_matches_by_query_pos(hit):
        if not is_match:
            yield qpos


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument(
        "type_strain_fasta", metavar="type-strain-fasta",
        type=argparse.FileType("r"),
        help="Type strain sequences FASTA file")
    p.add_argument(
        "reference_fasta", metavar="reference-fasta",
        help="Full-length reference sequences FASTA file")
    p.add_argument(
        "output_file", metavar="output-file", type=argparse.FileType("w"),
        help="Otuput file path")

    p.add_argument(
        "--batch-size",
        type=int, default=10,
        help=(
            "Number of query sequences to search simultaneously "
            "(default: %(default)s)"))
    p.add_argument(
        "--num-cpus", type=int,
        help="Number of CPUs to use in search (default: all the CPUs)")
    args = p.parse_args(argv)


    app = MismatchLocationApp(
        args.type_strain_fasta, args.reference_fasta, args.output_file,
        batch_size=args.batch_size, num_cpus=args.num_cpus,
    )
    app.run()

