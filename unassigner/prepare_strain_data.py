import argparse
import sys
import os

from unassigner.download import (
    get_url, clean,
    LTP_METADATA_URL, LTP_SEQS_URL,
    GG_SEQS_URL, GG_ACCESSIONS_URL,
    process_ltp_seqs, process_greengenes_seqs,
    )


def use_or_download(optional_fp, url, db_dir):
    if optional_fp is None:
        return get_url(url, os.path.join(db_dir, os.path.basename(url)))
    else:
        return optional_fp


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--ltp_metadata_fp", help=(
        "Filepath for LTP metadata (.csv file) "
        "[default:download from LTP website]"))
    p.add_argument("--ltp_seqs_fp", help=(
        "Filepath for unaligned 16S sequences from LTP (.fasta file) "
        "[default: download from LTP website]"))
    p.add_argument("--download_greengenes", action="store_true",
        help="Download GreenGenes reference files [default: False]")
    p.add_argument("--greengenes_accessions_fp", help=(
        "Filepath for table of GreenGenes accession numbers "
        "(.txt or .txt.gz file) "
        "[default: download from GreenGenes mirror]"))
    p.add_argument("--greengenes_seqs_fp", help=(
        "Filepath for unaligned 16S reference sequences "
        "(.fasta or .fasta.gz file) "
        "[default: download from GreenGenes mirror]"))
    p.add_argument("--clean", action="store_true", help=(
        "Remove all downloaded and processed files."))
    p.add_argument("--db-dir", help=(
        "Filepath to download the files to."))
    args = p.parse_args(argv)

    if args.db_dir:
        db_dir = args.db_dir
        if not os.path.exists(args.db_dir) and not args.clean:
            os.mkdir(args.db_dir)
    else:
        db_dir = os.getcwd()

    if args.clean is True:
        clean(db_dir)
        sys.exit(0)

    ltp_metadata_fp = use_or_download(
        args.ltp_metadata_fp, LTP_METADATA_URL, db_dir)
    ltp_seqs_fp = use_or_download(
        args.ltp_seqs_fp, LTP_SEQS_URL, db_dir)
    process_ltp_seqs(ltp_seqs_fp, db_dir)

    if args.download_greengenes:
        gg_seqs_fp = use_or_download(
            args.greengenes_seqs_fp, GG_SEQS_URL, db_dir)
        gg_accessions_fp = use_or_download(
            args.greengenes_accessions_fp, GG_ACCESSIONS_URL, db_dir)
        process_greengenes_seqs(gg_seqs_fp, gg_accessions_fp, db_dir)

def download_type_strain_data(output_dir=None, metadata_fp=None, seqs_fp=None):
    if output_dir is None:
        output_dir = os.getcwd()
    metadata_fp = use_or_download(metadata_fp, LTP_METADATA_URL, output_dir)
    seqs_fp = use_or_download(seqs_fp, LTP_SEQS_URL, output_dir)
    return process_ltp_seqs(seqs_fp, output_dir)
