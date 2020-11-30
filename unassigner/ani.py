import argparse
import collections
import os.path
import random
import re
import tempfile

from unassigner.download import get_url
from unassigner.align import VsearchAligner

class Refseq16SDatabase:
    fasta_fp = "refseq_16S.fasta"
    accession_fp = "refseq_16S_accessions.txt"
    hits_fp = "refseq_16S_hits.txt"

    def __init__(self):
        self.seqs = []
        self.seqid_to_assembly = {}

    def add_assembly(cls, assembly):
        assembly.download_rrna_16S()
        for seqid, seq in assembly.rrna_16S_seqs:
            self.seqid_to_assembly[seqid] = assembly
            self.seqs.append((seqid, seq))

    def save(self):
        with open(fasta_fp, "w") as f:
            write_fasta(f, self.seqs)
        with open(self.accession_fp, "w") as f:
            for seqid, assembly in self.seqid_to_assembly.items():
                f.write("{0}\t{1}\n".format(seqid, assembly.accession))

    def compute_pctids(self, min_pctid=97.0, threads=None):
        aligner = VsearchAligner(self.fasta_fp)
        hits = aligner.search(
            seqs=self.seqs, output_fp=self.hits_fp,
            min_id = min_pctid, threads=threads)
        for hit in hits:
            query = self.seqid_to_assembly[hit["qseqid"]]
            subject = self.seqid_to_assembly[hit["sseqid"]]
            pctid = hit["pctid"]
            yield AssemblyPair(query, subject, pctid)


class RefseqAssembly:
    summary_url = (
        "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/"
        "bacteria/assembly_summary.txt"
        )
    summary_fp = "refseq_bacteria_assembly_summary.txt"
    genome_dir = "genomes"
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
        self.seqs_16S = None
        for key, val in kwargs.items():
            setattr(self, key, val)

    @classmethod
    def parse_summary(cls, f):
        for line in f:
            line = line.strip()
            if line.startswith("#") or (line == ""):
                continue
            toks = line.split("\t")
            vals = dict(zip(cls.summary_cols, toks))
            accession = vals["assembly_accession"]
            yield cls(**vals)

    @classmethod
    def load(cls):
        if not os.path.exists(cls.assembly_summary_fp):
            get_url(cls.summary_url, cls.summary_fp)
        with open(cls.assembly_summary_fp, "r") as f:
            return {a.accession: a for a in cls.parse_summary(f)}

    def download_16S(self):
        if self.seqs_16S is not None:
            return
        with tempfile.NamedTemporaryFile() as f:
            get_url(self.rna_url, f.name)
            f.seek(0)
            seqs = list(parse_fasta(f))
        self.seqs_16S = [
            (desc, seq) for (desc, seq) in seqs if is_16S(desc)]

    @property
    def base_url(self):
        return re.sub("^ftp://", "https://", self.ftp_path)

    @property
    def rna_url(self):
        return "{0}/{1}_{2}_rna_from_genomic.fna.gz".format(
            self.base_url, self.accession, self.asm_name)

    @property
    def genome_fna_url(self):
        return "{0}/{1}_{2}_genomic.fna.gz".format(
            self.base_url, self.accession, self.asm_name)

    @property
    def genome_fna_fp(self):
        genome_filename = "{0}_{1}_genomic.fna.gz".format(
            self.accession, self.asm_name)
        return os.path.join(self.genome_dir, genome_filename)

    def download_genome_fna(self):
        if os.path.exists(self.genome_fna_fp):
            return
        if not os.path.exists(self.genome_dir):
            os.mkdir(self.genome_dir)
        get_url(self.genome_fna_url, self.genome_fna_fp)

def is_16S(desc):
    return "product=16S ribosomal RNA" in desc

def subsample_by(xs, fcn, n):
    groups = group_by(xs, fcn)
    return flatten(
        random.sample(group, n) for group in groups)


def flatten(xss):
    return [x for xs in xss for x in xs]


def group_by(xs, fcn):
    groups = collections.defaultdict(list)
    for x in xs:
        group = fcn(x)
        groups[group].append(x)
    return groups.values()


class AssemblyPair:
    def __init__(self, query, subject, pctid=None):
        self.query = query
        self.subject = subject
        self._pctid = pctid
        self.ani = None

    @property
    def pctid(self):
        return round(float(self._pctid), 1)

    def compute_ani(self):
        self.query.download_genome_fna()
        self.subject.download_genome_fna()
        # pyani get ANI for pair
        # save ANI to self.ani

    def format_output(self):
        return "{0}\t{1}\t{2}\t{3}\n".format(
            self.query.accession, self.subject.accession,
            self.pctid, self.ani,
        )


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument(
        "--output-file", type=argparse.FileType("w"),
        default="refseq_pctid_ani.tsv",
        help="Output file",
    )
    p.add_argument(
        "--min_pctid", type=float, default=97.0,
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

    # Stage 1: download RefSeq genome list
    assemblies = RefseqAssembly.load()

    # Stage 2: download 16S gene sequences
    pctid_db = Refseq16SDatabase()
    for assembly in assemblies:
        pctid_db.add_assembly(assembly)
    pctid_db.save()

    # Stage 3: compute 16S pctid
    all_pairs = rna_db.compute_pctids(args.min_pctid, args.threads)

    # Stage 4: select genome pairs for ANI
    random.seed(args.seed)
    selected_pairs = subsample_by(
        all_pairs, operator.attrgetter("pctid"), args.num_ani_pairs)

    # Stage 5: download full genomes for ANI comparison
    # Stage 6: compute ANI values
    # Stage 7: write final output file
    for pair in selected_pairs:
        pair.compute_ani()
        args.output_file.write(pair.format_output())
