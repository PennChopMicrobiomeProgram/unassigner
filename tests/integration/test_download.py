import pytest

from unassigner.download import LTP_VERSION, LTP_SEQS_URL, get_url
from unassigner.parse import parse_desc


@pytest.fixture
def download_LTP_fasta(tmp_path):
    fasta_fp = tmp_path / f"LTP_{LTP_VERSION}_blastdb.fasta"
    get_url(LTP_SEQS_URL, fasta_fp)
    yield fasta_fp


def test_parse_desc(download_LTP_fasta):
    with open(download_LTP_fasta) as f:
        line = f.readline()
        print(line)
        accession, species_name = parse_desc(line)

    assert accession
    assert species_name
