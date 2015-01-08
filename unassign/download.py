import os
import subprocess
import tarfile
import urllib2

from unassign.parse import parse_fasta, parse_greengenes_accessions

LTP_METADATA_COLS = [
    "accession",
    "start",
    "stop",
    "unknown",
    "fullname_ltp",
    "type_ltp",
    "hi_tax_ltp",
    "riskgroup_ltp",
    "url_lpsn_ltp",
    "tax_ltp",
    "rel_ltp",
    "NJ_support_pk4_ltp"
    ]
LTP_METADATA_URL = \
    "http://www.arb-silva.de/fileadmin/silva_databases/living_tree/LTP_release_119/LTPs119_SSU.csv"
LTP_SEQS_URL = \
    "http://www.arb-silva.de/fileadmin/silva_databases/living_tree/LTP_release_119/LTPs119_SSU.compressed.fasta"
GG_SEQS_URL = \
    "ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_5.fasta.gz"
GG_ACCESSIONS_URL = \
    "ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_5_accessions.txt.gz"
SPECIES_FASTA_FP = "species.fasta"
REFSEQS_FASTA_FP = "refseqs.fasta"


def clean():
    fps = [
        url_fp(LTP_METADATA_URL),
        url_fp(LTP_SEQS_URL),
        SPECIES_FASTA_FP,
        gunzip_fp(url_fp(GG_SEQS_URL)),
        gunzip_fp(url_fp(GG_ACCESSIONS_URL)),
        REFSEQS_FASTA_FP,
        ]
    fps += blastdb_fps(SPECIES_FASTA_FP)
    fps += blastdb_fps(REFSEQS_FASTA_FP)
    for fp in fps:
        if os.path.exists(fp):
            os.remove(fp)


def url_fp(url):
    return url.split('/')[-1]


def gunzip_fp(fp):
    return fp[:-3]


def blastdb_fps(fp):
    return [fp + ".nhr", fp + ".nin", fp + ".nsq"]


def get_url(url):
    subprocess.check_call(["wget", url])
    return url_fp(url)


def process_ltp_seqs(input_fp, output_fp=SPECIES_FASTA_FP):
    # Re-format FASTA file
    with open(input_fp) as f_in:
        seqs = parse_fasta(f_in)
        with open(output_fp, "w") as f_out:
            for desc, seq in seqs:
                name = desc.split("\t")[0]
                f_out.write(">%s\n%s\n" % (name, seq))

    # Make BLAST database
    subprocess.check_call([
        "makeblastdb",
        "-dbtype", "nucl",
        "-in", output_fp,
        ])

    return output_fp


def process_greengenes_seqs(seqs_fp, accessions_fp, output_fp=REFSEQS_FASTA_FP):
    # Extract table of accessions
    subprocess.check_call(["gunzip", "-f", accessions_fp])
    accessions_fp = gunzip_fp(accessions_fp)

    # Load accessions
    gg_accessions = {}
    with open(accessions_fp) as f:
        for ggid, src, acc in parse_greengenes_accessions(f):
            gg_accessions[ggid] = (acc, src)

    # Extract FASTA file
    subprocess.check_call(["gunzip", "-f", seqs_fp])
    seqs_fp = gunzip_fp(seqs_fp)

    # Re-label seqs with accession numbers
    with open(seqs_fp) as f_in:
        seqs = parse_fasta(f_in)
        with open(output_fp, "w") as f_out:
            for ggid, seq in seqs:
                acc, src = gg_accessions[ggid]
                f_out.write(">%s %s %s\n%s\n" % (acc, src, ggid, seq))
    
    # Make BLAST database
    subprocess.check_call([
        "makeblastdb",
        "-dbtype", "nucl",
        "-in", output_fp,
        ])

    return output_fp
