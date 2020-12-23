import argparse
import collections
import logging
import os
import os.path
import random
import re
import shutil
import subprocess
import urllib.error

from unassigner.download import get_url
from unassigner.parse import parse_fasta, write_fasta
from unassigner.align import DEFAULT_BLAST_FIELDS, BLAST_TO_VSEARCH, VsearchAligner

class Refseq16SDatabase:
    def __init__(self, working_dir=None):
        if working_dir is None:
            working_dir = os.getcwd()

        self.working_dir = working_dir
        self.ssu_accession_fp = os.path.join(
            self.working_dir, "refseq_16S_accessions_all.txt")
        self.ssu_fasta_fp = os.path.join(
            self.working_dir, "refseq_16S_all.fasta")
        self.query_fp = os.path.join(
            self.working_dir, "temp_query.fasta")
        self.hits_fp = os.path.join(
            self.working_dir, "temp_hits.txt")

        self.seqs = {}
        self.assemblies = {}
        self.seqids_by_assembly = collections.defaultdict(list)

    def load(self, assemblies):
        was_loaded_from_file = self._load_from_files(assemblies)
        if not was_loaded_from_file:
            for assembly in assemblies:
                self._add_assembly(assembly)
            self._save_to_files()
        return was_loaded_from_file

    def _load_from_files(self, assemblies):
        if os.path.exists(self.ssu_accession_fp):
            with open(self.ssu_accession_fp) as f:
                self._load_accessions(f, assemblies)
            with open(self.ssu_fasta_fp) as f:
                self._load_seqs(f)
            return True
        return False

    def _save_to_files(self):
        with open(self.ssu_accession_fp, "w") as f:
            self._save_accessions(f)
        with open(self.ssu_fasta_fp, "w") as f:
            self._save_seqs(f)

    def _add_assembly(self, assembly):
        seqs = list(assembly.ssu_seqs)
        # Avoid writing duplicate genes for the same genome
        seen = set()
        for desc, seq in seqs:
            if seq not in seen:
                logging.info(desc)
                seqid = desc.split()[0]
                self.assemblies[seqid] = assembly
                self.seqs[seqid] = seq
                self.seqids_by_assembly[assembly.accession].append(seqid)
                seen.add(seq)

    def _save_accessions(self, f):
        for seqid, assembly in self.assemblies.items():
            f.write("{0}\t{1}\n".format(seqid, assembly.accession))

    def _load_accessions(self, f, assemblies):
        assemblies_by_accession = {a.accession: a for a in assemblies}
        for line in f:
            toks = line.strip().split()
            seqid = toks[0]
            accession = toks[1]
            assembly = assemblies_by_accession[accession]
            self.assemblies[seqid] = assembly
            self.seqids_by_assembly[accession].append(seqid)

    def _save_seqs(self, f):
        write_fasta(f, self.seqs.items())

    def _load_seqs(self, f):
        for seqid, seq in parse_fasta(f):
            self.seqs[seqid] = seq

    def search_at_pctid(self, query_seqid, pctid, threads=None):
        pctid_str = "{:.1f}".format(pctid)
        print("Searching {0} at {1} pct identity".format(
            query_seqid, pctid_str))

        # Keep previous round of files for debugging
        if os.path.exists(self.query_fp):
            os.rename(self.query_fp, self.query_fp + "_prev")
        if os.path.exists(self.hits_fp):
            os.rename(self.hits_fp, self.hits_fp + "_prev")

        query_seq = self.seqs[query_seqid]
        query_seqs = [(query_seqid, query_seq)]

        # Must set minimum id a bit lower for the search
        min_pctid = pctid - 0.1
        min_id = min_pctid / 100
        fields = ["qseqid", "sseqid", "pident"]

        search_params = {
            "min_id": min_id,
            "maxaccepts": 10000,
            "threads": threads,
        }
        vsearch_aligner = VsearchAligner(self.ssu_fasta_fp)
        vsearch_aligner.fields = ["qseqid", "sseqid", "pident"]
        vsearch_aligner.convert_types = False
        hits = vsearch_aligner.search(
            query_seqs, self.query_fp, self.hits_fp, search_params)
        for hit in hits:
            if hit["qseqid"] == hit["sseqid"]:
                continue
            if hit["pident"] == pctid_str:
                query = self.assemblies[hit["qseqid"]]
                subject = self.assemblies[hit["sseqid"]]
                pctid = hit["pident"]
                yield AssemblyPair(
                    query, subject, pctid,
                    hit["qseqid"], hit["sseqid"])


class PctidAligner:
    field_names = ["qseqid", "sseqid", "pident"]

    def __init__(self, ref_seqs_fp):
        self.ref_seqs_fp = ref_seqs_fp

    @property
    def ref_seqs_udb_fp(self):
        base_fp, _ = os.path.splitext(self.ref_seqs_fp)
        return base_fp + ".udb"

    def make_reference_udb(self):
        if os.path.exists(self.ref_seqs_udb_fp):
            return None
        args = [
            "vsearch",
            "--makeudb_usearch", self.ref_seqs_fp,
            "--output", self.ref_seqs_udb_fp,
        ]
        return subprocess.check_call(args)

    def search(
            self, query_fp, output_fp,
            min_id=0.970, fields=DEFAULT_BLAST_FIELDS,
            threads=None, maxaccepts=10000):
        self.make_reference_udb()
        min_id_arg = "{:.3f}".format(min_id)
        vsearch_fields = [BLAST_TO_VSEARCH[f] for f in fields]
        vsearch_fields_arg = "+".join(vsearch_fields)
        maxaccepts_arg = "{:d}".format(maxaccepts)
        args =[
            "vsearch", "--usearch_global", query_fp,
            "--db", self.ref_seqs_udb_fp,
            "--iddef", "2",
            "--id", min_id_arg,
            "--userout", output_fp,
            "--userfields", vsearch_fields_arg,
            "--maxaccepts", maxaccepts_arg,
            ]
        if threads is not None:
            threads_arg = "{:d}".format(threads)
            args.extend(["--threads", threads_arg])
        subprocess.check_call(args, stderr=subprocess.DEVNULL)

    def parse(self, f):
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            vals = line.split("\t")
            hit = dict(zip(self.field_names, vals))
            if hit["qseqid"] != hit["sseqid"]:
                yield hit


class RefseqAssembly:
    summary_url = (
        "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/"
        "bacteria/assembly_summary.txt"
        )
    summary_cols = [
        "assembly_accession", "bioproject", "biosample", "wgs_master",
        "refseq_category", "taxid", "species_taxid", "organism_name",
        "infraspecific_name", "isolate", "version_status",
        "assembly_level", "release_type", "genome_rep", "seq_rel_date",
        "asm_name", "submitter", "gbrs_paired_asm", "paired_asm_comp",
        "ftp_path", "excluded_from_refseq", "relation_to_type_material"
    ]

    def __init__(self, assembly_accession, ftp_path, **kwargs):
        self.accession = assembly_accession
        self.ftp_path = ftp_path
        for key, val in kwargs.items():
            setattr(self, key, val)
        self._ssu_seqs = None
        self.genome_dir = "genome_fasta"
        self.rna_dir = "rna_fasta"

    @property
    def _base_url(self):
        return re.sub("^ftp://", "https://", self.ftp_path)

    @property
    def _basename(self):
        return os.path.basename(self.ftp_path)

    @classmethod
    def parse_summary(cls, f):
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#") or (line == ""):
                continue
            toks = line.split("\t")
            vals = dict(zip(cls.summary_cols, toks))
            if vals["ftp_path"] == "na":
                continue
            yield cls(**vals)

    @property
    def ssu_seqs(self):
        if self._ssu_seqs is None:
            self._download_rna()
            self._load_ssu_seqs()
        return self._ssu_seqs

    @property
    def rna_url(self):
        return "{0}/{1}_rna_from_genomic.fna.gz".format(
            self._base_url, self._basename)

    @property
    def rna_fp(self):
        rna_filename = "{0}_rna_from_genomic.fna".format(self._basename)
        return os.path.join(self.rna_dir, rna_filename)

    def _download_rna(self):
        if os.path.exists(self.rna_fp):
            return
        if not os.path.exists(self.rna_dir):
            os.mkdir(self.rna_dir)
        try:
            logging.info(
                "Downloading RNA file for {0}".format(self.accession))
            get_url(self.rna_url, self.rna_fp + ".gz")
            subprocess.check_call(["gunzip", "-q", self.rna_fp + ".gz"])
        except Exception as e:
            logging.warning(
                "Error downloading RNA file for {0}:\n{1}".format(
                    self.accession, e))

    def _load_ssu_seqs(self):
        if os.path.exists(self.rna_fp):
            with open(self.rna_fp, "rt") as f:
                seqs = list(parse_fasta(f))
            res = [(desc, seq) for (desc, seq) in seqs if is_16S(desc)]
            self._ssu_seqs = res

    @property
    def genome_url(self):
        return "{0}/{1}_genomic.fna.gz".format(
            self._base_url, self._basename)

    @property
    def genome_fp(self):
        genome_filename = "{0}_genomic.fna.gz".format(self._basename)
        return os.path.join(self.genome_dir, genome_filename)

    def download_genome(self, genome_dir=None, genome_fname=None):
        if genome_dir is None:
            genome_dir = self.genome_dir
        if genome_fname is None:
            genome_fp = self.genome_fp
        else:
            genome_fp = os.path.join(genome_dir, genome_fname)
        if os.path.exists(genome_fp):
            return
        if not os.path.exists(genome_dir):
            os.mkdir(genome_dir)
        get_url(self.genome_url, genome_fp)


def is_16S(desc):
    return "product=16S ribosomal RNA" in desc


def remove_files(target_dir):
    if os.path.exists(target_dir):
        for filename in os.listdir(target_dir):
            os.remove(os.path.join(target_dir, filename))


class AssemblyPair:
    genome_dir = "temp_genomes"
    ani_dir = "temp_ani"
    cache_dir = "temp_genomes_cache"
    ani_cache = {}

    def __init__(
            self, query, subject, pctid=None,
            query_seqid=None, subject_seqid=None):
        self.query = query
        self.subject = subject
        self._pctid = pctid
        self.ani = None
        self.query_seqid = query_seqid
        self.subject_seqid = subject_seqid

    @property
    def pctid(self):
        return round(float(self._pctid), 1)

    def compute_ani(self):
        ani_key = (self.query.accession, self.subject.accession)
        if ani_key in self.ani_cache:
            print("ANI is cached")
            self.ani = self.ani_cache[ani_key]
            return

        remove_files(self.genome_dir)
        
        query_fname = "{0}.fna.gz".format(self.query.accession)
        query_genome_fp = os.path.join(self.genome_dir, query_fname)
        query_cache_fp = os.path.join(self.cache_dir, query_fname)
        if os.path.exists(query_cache_fp):
            os.rename(query_cache_fp, query_genome_fp)
        else:    
            self.query.download_genome(self.genome_dir, query_fname)
            shutil.copy(query_genome_fp, query_cache_fp)
        subprocess.check_call(["gunzip", "-f", query_genome_fp])

        subject_fname = "{0}.fna.gz".format(self.subject.accession)
        subject_genome_fp = os.path.join(self.genome_dir, subject_fname)
        self.subject.download_genome(self.genome_dir, subject_fname)
        subprocess.check_call(["gunzip", "-f", subject_genome_fp])
        
        # pyani get ANI for pair
        shutil.rmtree(self.ani_dir)
        # save ANI to self.ani
        print(
            "Computing ANI for", self.query.accession, "and",
            self.subject.accession)
        subprocess.check_call([
            "average_nucleotide_identity.py",
            "-i", self.genome_dir,
            "-o", self.ani_dir,
        ])

        ani_fp = os.path.join(self.ani_dir, "ANIm_percentage_identity.tab")
        with open(ani_fp) as f:
            ani = parse_pairwise_ani(f)
        print("ANI:", ani)
        self.ani_cache[ani_key] = ani
        self.ani = ani

    def format_output(self):
        if self.query_seqid and self.subject_seqid:
            return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                self.query.accession, self.subject.accession,
                self.query_seqid, self.subject_seqid,
                self.pctid, self.ani,
            )
        return "{0}\t{1}\t{2}\t{3}\n".format(
            self.query.accession, self.subject.accession,
            self.pctid, self.ani,
        )

def parse_pairwise_ani(f):
    # first line is a header
    next(f)
    # second line has ani as the third token
    second_line = next(f)
    toks = second_line.strip().split()
    return toks[2]

def pctid_range(min_pctid):
    assert 50.0 < min_pctid <= 100.0
    current_val = 100.0
    while current_val > min_pctid:
        yield current_val
        current_val = current_val - 0.1

def main_sampling(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument(
        "--output-file", type=argparse.FileType("w"),
        default="refseq_pctid_ani.tsv",
        help="Output file",
    )
    p.add_argument(
        "--assembly-summary-file",
        default="refseq_bacteria_assembly_summary.txt",
        help=(
            "File to use for RefSeq assembly summary. Will be downloaded if "
            "not already present."),
    )
    p.add_argument(
        "--min-pctid", type=float, default=97.0,
        help="Minimum 16S percent ID",
    )
    p.add_argument(
        "--num-threads", type=int,
        help="Number of threads for 16S percent ID (default: use all CPUs)",
    )
    p.add_argument(
        "--num-ani", type=int, default=100,
        help="Number of genome pairs on which to evaluate ANI",
    )
    p.add_argument(
        "--seed", type=int, default=42,
        help="Random number seed",
    )
    args = p.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(message)s",
    )

    args.output_file.write(
        "query_assembly\tsubject_assembly\t"
        "query_seqid\tsubject_seqid\t"
        "pctid\tani\n")

    random.seed(args.seed)

    if not os.path.exists(args.assembly_summary_file):
        get_url(RefseqAssembly.summary_url, args.assembly_summary_file)

    with open(args.assembly_summary_file, "r") as f:
        assemblies = list(RefseqAssembly.parse_summary(f))

    db = Refseq16SDatabase()
    db.load(assemblies)

    pctid_vals = list(pctid_range(args.min_pctid)) * args.num_ani

    for current_pctid in pctid_vals:
        found = False
        while not found:
            query_assembly = random.choice(assemblies)
            query_assembly_seqids = db.seqids_by_assembly[query_assembly.accession]
            if not query_assembly_seqids:
                continue
            query_seqid = random.choice(query_assembly_seqids)
            assembly_pairs = db.search_at_pctid(
                query_seqid, current_pctid, threads=args.num_threads)
            assembly_pairs = list(assembly_pairs)
            if assembly_pairs:
                selected_pair = random.choice(assembly_pairs)
                try:
                    selected_pair.compute_ani()
                    args.output_file.write(selected_pair.format_output())
                    args.output_file.flush()
                except Exception as e:
                    logging.warning(
                        "Exception during ANI computation:\n{0}".format(e))
                else:
                    found = True

