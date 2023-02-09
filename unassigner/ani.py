import argparse
import collections
import os
import os.path
import random
import re
import shutil
import subprocess
import urllib.error

from unassigner.download import get_url
from unassigner.parse import parse_fasta, write_fasta


class Refseq16SDatabase:
    def __init__(
        self, fasta_fp="refseq_16S.fasta", accession_fp="refseq_16S_accessions.txt"
    ):
        self.fasta_fp = fasta_fp
        self.accession_fp = accession_fp
        self.seqs = {}
        self.assemblies = {}
        self.seqids_by_assembly = collections.defaultdict(list)

    def add_assembly(self, assembly, select_random=False):
        seqs = list(assembly.ssu_seqs)
        if select_random and (len(seqs) > 0):
            seqs = [random.choice(seqs)]

        # Avoid writing duplicate genes for the same genome
        seen = set()
        for desc, seq in seqs:
            if seq not in seen:
                print(desc)
                seqid = desc.split()[0]
                self.assemblies[seqid] = assembly
                self.seqs[seqid] = seq
                self.seqids_by_assembly[assembly.accession].append(seqid)
                seen.add(seq)

    def load(self, assemblies):
        with open(self.accession_fp, "r") as f:
            for line in f:
                toks = line.strip().split()
                seqid = toks[0]
                accession = toks[1]
                assembly = assemblies[accession]
                self.assemblies[seqid] = assembly
                self.seqids_by_assembly[assembly.accession].append(seqid)
        with open(self.fasta_fp, "r") as f:
            for seqid, seq in parse_fasta(f):
                self.seqs[seqid] = seq

    def save(self):
        with open(self.fasta_fp, "w") as f:
            write_fasta(f, self.seqs.items())
        with open(self.accession_fp, "w") as f:
            for seqid, assembly in self.assemblies.items():
                f.write("{0}\t{1}\n".format(seqid, assembly.accession))

    def compute_pctids(self, min_pctid=97.0, threads=None):
        aligner = PctidAligner(self.fasta_fp)
        if not os.path.exists(aligner.hits_fp):
            aligner.search(min_pctid=min_pctid, threads=threads)
        with open(aligner.hits_fp) as f:
            for hit in aligner.parse(f):
                query = self.assemblies[hit["qseqid"]]
                subject = self.assemblies[hit["sseqid"]]
                pctid = hit["pident"]
                yield AssemblyPair(query, subject, pctid)

    def search_one(self, query_seqid, pctid, threads=None):
        pctid_str = "{:.1f}".format(pctid)
        print("Searching", query_seqid, "at", pctid_str, "pct identity")
        query_seq = self.seqs[query_seqid]
        query_fp = "temp_query.fasta"
        if os.path.exists(query_fp):
            os.rename(query_fp, "temp_prev_query.fasta")
        query_hits_fp = "temp_query_hits.txt"
        if os.path.exists(query_hits_fp):
            os.rename(query_hits_fp, "temp_prev_query_hits.txt")
        with open(query_fp, "w") as f:
            write_fasta(f, [(query_seqid, query_seq)])
        aligner = PctidAligner(self.fasta_fp)
        aligner.search(
            query_fp, query_hits_fp, min_pctid=pctid, threads=threads, max_hits=10000
        )
        with open(query_hits_fp) as f:
            hits = aligner.parse(f)
            for hit in hits:
                if hit["pident"] == pctid_str:
                    query = self.assemblies[hit["qseqid"]]
                    subject = self.assemblies[hit["sseqid"]]
                    pctid = hit["pident"]
                    yield AssemblyPair(
                        query, subject, pctid, hit["qseqid"], hit["sseqid"]
                    )

    def search_seq(self, query_seqid, query_seq, min_pctid=90.0, threads=None):
        query_fp = "temp_query.fasta"
        if os.path.exists(query_fp):
            os.rename(query_fp, "temp_prev_query.fasta")
        query_hits_fp = "temp_query_hits.txt"
        if os.path.exists(query_hits_fp):
            os.rename(query_hits_fp, "temp_prev_query_hits.txt")
        with open(query_fp, "w") as f:
            write_fasta(f, [(query_seqid, query_seq)])
        aligner = PctidAligner(self.fasta_fp)
        aligner.search(
            query_fp,
            query_hits_fp,
            min_pctid=min_pctid,
            threads=threads,
            max_hits=10000,
        )
        with open(query_hits_fp) as f:
            hits = aligner.parse(f)
            for hit in hits:
                query = self.assemblies[hit["qseqid"]]
                subject = self.assemblies[hit["sseqid"]]
                pctid = hit["pident"]
                if query.accession != subject.accession:
                    yield AssemblyPair(
                        query, subject, pctid, hit["qseqid"], hit["sseqid"]
                    )


class PctidAligner:
    field_names = ["qseqid", "sseqid", "pident"]
    hits_fp = "refseq_16S_hits.txt"

    def __init__(self, fasta_fp):
        self.fasta_fp = fasta_fp

    @property
    def reference_udb_fp(self):
        base_fp, _ = os.path.splitext(self.fasta_fp)
        return base_fp + ".udb"

    def make_reference_udb(self):
        if os.path.exists(self.reference_udb_fp):
            return None
        args = [
            "vsearch",
            "--makeudb_usearch",
            self.fasta_fp,
            "--output",
            self.reference_udb_fp,
        ]
        return subprocess.check_call(args)

    def search(
        self, input_fp=None, hits_fp=None, min_pctid=97.0, threads=None, max_hits=10000
    ):
        if input_fp is None:
            input_fp = self.fasta_fp
        if hits_fp is None:
            hits_fp = self.hits_fp
        self.make_reference_udb()
        # 97.0 --> 0.97
        min_id = "{:.3f}".format(min_pctid / 100)
        args = [
            "vsearch",
            "--usearch_global",
            input_fp,
            "--db",
            self.reference_udb_fp,
            "--userout",
            hits_fp,
            "--iddef",
            "2",
            "--id",
            min_id,
            "--userfields",
            "query+target+id2",
        ]
        if max_hits is None:
            args.extend(
                [
                    "--maxaccepts",
                    "0",
                    "--maxrejects",
                    "0",
                ]
            )
        else:
            args.extend(
                [
                    "--maxaccepts",
                    str(max_hits),
                ]
            )
        if threads is not None:
            args.extend(["--threads", str(threads)])
        print(args)
        subprocess.check_call(args)
        return hits_fp

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
        "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"
    )
    summary_fp = "refseq_bacteria_assembly_summary.txt"
    genome_dir = "genome_fasta"
    rna_dir = "rna_fasta"
    summary_cols = [
        "assembly_accession",
        "bioproject",
        "biosample",
        "wgs_master",
        "refseq_category",
        "taxid",
        "species_taxid",
        "organism_name",
        "infraspecific_name",
        "isolate",
        "version_status",
        "assembly_level",
        "release_type",
        "genome_rep",
        "seq_rel_date",
        "asm_name",
        "submitter",
        "gbrs_paired_asm",
        "paired_asm_comp",
        "ftp_path",
        "excluded_from_refseq",
        "relation_to_type_material",
    ]

    def __init__(self, assembly_accession, ftp_path, **kwargs):
        self.accession = assembly_accession
        self.ftp_path = ftp_path
        for key, val in kwargs.items():
            setattr(self, key, val)
        self._ssu_seqs = None

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
        if self._ssu_seqs is not None:
            return self._ssu_seqs
        if not os.path.exists(self.rna_fp):
            try:
                self.download_rna()
            except urllib.error.HTTPError as e:
                print(self.accession)
                print(e)
                return []
        with open(self.rna_fp, "rt") as f:
            seqs = list(parse_fasta(f))
        res = [(desc, seq) for (desc, seq) in seqs if is_16S(desc)]
        self._ssu_seqs = res
        return res

    @classmethod
    def load(cls):
        if not os.path.exists(cls.summary_fp):
            get_url(cls.summary_url, cls.summary_fp)
        with open(cls.summary_fp, "r") as f:
            return {a.accession: a for a in cls.parse_summary(f)}

    @property
    def base_url(self):
        return re.sub("^ftp://", "https://", self.ftp_path)

    @property
    def basename(self):
        return os.path.basename(self.ftp_path)

    @property
    def rna_url(self):
        return "{0}/{1}_rna_from_genomic.fna.gz".format(self.base_url, self.basename)

    @property
    def rna_fp(self):
        rna_filename = "{0}_rna_from_genomic.fna".format(self.basename)
        return os.path.join(self.rna_dir, rna_filename)

    def download_rna(self):
        if os.path.exists(self.rna_fp):
            return
        if not os.path.exists(self.rna_dir):
            os.mkdir(self.rna_dir)
        print("Downloading 16S seqs for ", self.accession)
        get_url(self.rna_url, self.rna_fp + ".gz")
        subprocess.check_call(["gunzip", "-q", self.rna_fp + ".gz"])

    @property
    def genome_url(self):
        return "{0}/{1}_genomic.fna.gz".format(self.base_url, self.basename)

    @property
    def genome_fp(self):
        genome_filename = "{0}_genomic.fna.gz".format(self.basename)
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


def subsample_by(xs, fcn, n):
    groups = group_by(xs, fcn)
    for group in groups:
        if len(group) > n:
            group = random.sample(group, n)
        for x in group:
            yield x


def flatten(xss):
    return [x for xs in xss for x in xs]


def group_by(xs, fcn):
    groups = collections.defaultdict(list)
    for x in xs:
        group = fcn(x)
        groups[group].append(x)
    return groups.values()


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
        self, query, subject, pctid=None, query_seqid=None, subject_seqid=None
    ):
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
        print("Computing ANI for", self.query.accession, "and", self.subject.accession)
        subprocess.check_call(
            [
                "average_nucleotide_identity.py",
                "-i",
                self.genome_dir,
                "-o",
                self.ani_dir,
            ]
        )

        ani_fp = os.path.join(self.ani_dir, "ANIm_percentage_identity.tab")
        with open(ani_fp) as f:
            ani = parse_pairwise_ani(f)
        print("ANI:", ani)
        self.ani_cache[ani_key] = ani
        self.ani = ani

    def format_output(self):
        if self.query_seqid and self.subject_seqid:
            return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                self.query.accession,
                self.subject.accession,
                self.query_seqid,
                self.subject_seqid,
                self.pctid,
                self.ani,
            )
        return "{0}\t{1}\t{2}\t{3}\n".format(
            self.query.accession,
            self.subject.accession,
            self.pctid,
            self.ani,
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
        "--output-file",
        type=argparse.FileType("w"),
        default="refseq_pctid_ani.tsv",
        help="Output file",
    )
    p.add_argument(
        "--min_pctid",
        type=float,
        default=97.0,
        help="Minimum 16S percent ID",
    )
    p.add_argument(
        "--num-threads",
        type=int,
        help="Number of threads for 16S percent ID (default: use all CPUs)",
    )
    p.add_argument(
        "--num-ani",
        type=int,
        default=100,
        help="Number of genome pairs on which to evaluate ANI",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random number seed",
    )
    args = p.parse_args()
    args.output_file.write(
        "query_assembly\tsubject_assembly\t"
        "query_seqid\tsubject_seqid\t"
        "pctid\tani\n"
    )

    # Set seed for 16S selection
    random.seed(args.seed)

    # Load all assemblies
    assemblies = RefseqAssembly.load()

    # 16S database with one random sequence from each assembly
    db = Refseq16SDatabase("refseq_16S_all.fasta", "refseq_16S_accessions_all.txt")
    if os.path.exists(db.accession_fp):
        db.load(assemblies)
    else:
        for assembly in assemblies.values():
            db.add_assembly(assembly)
        db.save()

    pctid_vals = list(pctid_range(args.min_pctid)) * args.num_ani

    # Set seed again
    random.seed(args.seed + 1)
    assembly_list = list(assemblies.values())
    for current_pctid in pctid_vals:
        found = False
        while not found:
            # randomly select assembly
            query_assembly = random.choice(assembly_list)
            # randomly select query sequence from assembly
            query_assembly_seqids = db.seqids_by_assembly[query_assembly.accession]
            # next loop if we don't have any sequences for this assembly
            if not query_assembly_seqids:
                continue
            query_seqid = random.choice(query_assembly_seqids)
            assembly_pairs = db.search_one(
                query_seqid, current_pctid, threads=args.num_threads
            )
            assembly_pairs = list(assembly_pairs)
            if assembly_pairs:
                try:
                    # randomly select one result
                    selected_pair = random.choice(assembly_pairs)
                    selected_pair.compute_ani()
                    args.output_file.write(selected_pair.format_output())
                    args.output_file.flush()
                except Exception as e:
                    print(e)
                else:
                    found = True
