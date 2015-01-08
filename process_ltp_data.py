#!/usr/bin/env python
import optparse
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
    p = optparse.OptionParser()
    p.add_option("--ltp_metadata_fp", help=(
        "Filepath for LTP metadata (.csv file) "
        "[default:download from LTP website]"))
    p.add_option("--ltp_seqs_fp", help=(
        "Filepath for unaligned 16S sequences from LTP (.fasta file) "
        "[default: download from LTP website]"))
    p.add_option("--greengenes_accessions_fp", help=(
        "Filepath for table of GreenGenes accession numbers (.txt.gz file) "
        "[default: download from GreenGenes mirror]"))
    p.add_option("--greengenes_seqs_fp", help=(
        "Filepath for unaligned 16S reference sequences (.fasta.gz file) "
        "[default: download from GreenGenes mirror]"))
    p.add_option("--clean", action="store_true", help=(
        "Remove all downloaded and processed files."))
    opts, args = p.parse_args(argv)

    if opts.clean is True:
        clean()
        sys.exit(0)

    ltp_metadata_fp = use_or_download(
        opts.ltp_metadata_fp, LTP_METADATA_URL)
    ltp_seqs_fp = use_or_download(
        opts.ltp_seqs_fp, LTP_SEQS_URL)
    process_ltp_seqs(ltp_seqs_fp)

    gg_seqs_fp = use_or_download(
        opts.greengenes_seqs_fp, GG_SEQS_URL)
    gg_accessions_fp = use_or_download(
        opts.greengenes_accessions_fp, GG_ACCESSIONS_URL)
    process_greengenes_seqs(gg_seqs_fp, gg_accessions_fp)


if __name__ == "__main__":
    main()
