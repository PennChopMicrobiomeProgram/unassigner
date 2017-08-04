import argparse
import sys

from unassign.download import (
    get_url, clean,
    LTP_METADATA_URL, LTP_SEQS_URL,
    GG_SEQS_URL, GG_ACCESSIONS_URL,
    process_ltp_seqs, process_greengenes_seqs,
    )


def use_or_download(optional_fp, url):
    if optional_fp is None:
        return get_url(url)
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
    args = p.parse_args(argv)

    if args.clean is True:
        clean()
        sys.exit(0)

    ltp_metadata_fp = use_or_download(
        args.ltp_metadata_fp, LTP_METADATA_URL)
    ltp_seqs_fp = use_or_download(
        args.ltp_seqs_fp, LTP_SEQS_URL)
    process_ltp_seqs(ltp_seqs_fp)

    gg_seqs_fp = use_or_download(
        args.greengenes_seqs_fp, GG_SEQS_URL)
    gg_accessions_fp = use_or_download(
        args.greengenes_accessions_fp, GG_ACCESSIONS_URL)
    process_greengenes_seqs(gg_seqs_fp, gg_accessions_fp)
