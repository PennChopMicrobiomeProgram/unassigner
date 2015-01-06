#!/usr/bin/env python
import csv
import itertools
import optparse
import urllib2
import tarfile
import subprocess
from cStringIO import StringIO

LTP_BASE_URL = "http://www.arb-silva.de/fileadmin/silva_databases/living_tree/"
LTP_METADATA_URL = LTP_BASE_URL + "LTP_release_115/LTPs115_SSU.csv"
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
LTP_SEQS_URL = LTP_BASE_URL + "LTP_release_115/LTPs115.compressed.fasta.tar.gz"
GG_SEQS_URL = "ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_5.fasta.gz"

def process_metadata(f):
    rows = csv.reader(f, delimiter=";")
    for row in itertools.islice(rows, 0, 15):
        rec = dict(zip(LTP_METADATA_COLS, row))
        print rec

def parse_ltp_fasta(f):
    f = iter(f)
    desc = f.next().strip()[1:]
    seq = StringIO()
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            yield desc, seq.getvalue()
            desc = line[1:]
            seq = StringIO()
        else:
            seq.write(line.replace(" ", "").replace("U", "T"))
    yield desc, seq.getvalue()


BLAST_FMT = (
    "7 qseqid sseqid pident length mismatch gapopen "
    "qstart qend sstart send qlen slen qseq sseq")
    

def get_url(url):
    fp = url.split('/')[-1]
    local_file = open(fp, "wb")
    print "Downloading", url
    web_file = urllib2.urlopen(url)
    local_file.write(web_file.read())
    local_file.close()
    print "Done"
    return fp


def main(argv=None):
    p = optparse.OptionParser()
    p.add_option("--metadata_fp", help=(
        "Filepath for LTP metadata (.csv file) "
        "[default:download from LTP website]"))
    p.add_option("--type_strain_seqs_fp", help=(
        "Filepath for unaligned 16S sequences from LTP (.tar.gz file) "
        "[default: download from LTP website]"))
    p.add_option("--all_seqs_fp", help=(
        "Filepath for unaligned 16S reference sequences (.fasta.gz file) "
        "[default: download from GreenGenes mirror]"))
    
    p.add_option("--output_dir", help=(
        "Output directory"))
    opts, args = p.parse_args(argv)

    if opts.metadata_fp is None:
        metadata_file = urllib2.urlopen(LTP_METADATA_URL)
    else:
        metadata_file = open(opts.metadata_fp)
    metadata_output_file = open("species.csv", "w")
    metadata_output_file.write(metadata_file.read())
    metadata_output_file.close()
    #process_metadata(metadata_file)

    if opts.type_strain_seqs_fp is None:
        seqs_fp = get_url(LTP_SEQS_URL)
    else:
        seqs_fp = opts.type_strain_seqs_fp
    seqs_archive = tarfile.open(seqs_fp, mode="r:gz")
    seqs_member = seqs_archive.getnames()[0]
    seqs_file = seqs_archive.extractfile(seqs_member)
    seqs = parse_ltp_fasta(seqs_file)

    species_fp = "species.fasta"
    species_file = open(species_fp, "w")
    for desc, seq in seqs:
        name = desc.split("\t")[0]
        species_file.write(">%s\n%s\n" % (name, seq))
    species_file.close()
    subprocess.check_call([
        "makeblastdb",
        "-dbtype", "nucl",
        "-in", species_fp,
        ])

    if opts.all_seqs_fp is None:
        all_seqs_gz_fp = get_url(GG_SEQS_URL)
    else:
        all_seqs_gz_fp = opts.all_seqs_fp
    subprocess.check_call(["gunzip", all_seqs_gz_fp])
    all_seqs_fp = all_seqs_gz_fp[:-3]
    subprocess.check_call([
        "makeblastdb",
        "-dbtype", "nucl",
        "-in", all_seqs_fp,
        ])

    #species_blast_fp = "species_blast.txt"
    #args = [
    #    "blastn",
    #    "-evalue", "1e-5",
    #    "-max_target_seqs", "100",
    #    "-outfmt", BLAST_FMT,
    #    "-db", all_seqs_fp,
    #    "-query", species_fp,
    #    "-out", species_blast_fp,
    #    ]
    #subprocess.check_call(args)

    
if __name__ == "__main__":
    main()
