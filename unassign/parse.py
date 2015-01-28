from cStringIO import StringIO

def parse_fasta(f):
    """Parse a FASTA format file.

    Parameters
    ----------
    f : File object or iterator returning lines in FASTA format.

    Returns
    -------
    An iterator of tuples containing two strings
        First string is the sequence description, second is the
        sequence.

    Notes
    -----
    This function removes whitespace in the sequence and translates
    "U" to "T", in order to accommodate FASTA files downloaded from
    SILVA and the Living Tree Project.
    """
    f = iter(f)
    desc = next(f).strip()[1:]
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


def trim_seq_desc(seqs):
    for desc, seq in seqs:
        seqid = desc.split()[0]
        yield seqid, seq


def write_fasta(f, seqs):
    for desc, seq in seqs:
        f.write(">%s\n%s\n" % (desc, seq))


def load_fasta(filepath):
    """Load all sequences from a FASTA file

    Parameters
    ----------
    fasta_fp : Input filepath, FASTA format.

    Returns
    -------
    A dictionary mapping sequence identifiers to sequences.
    """
    with open(filepath) as f:
        seqs = trim_seq_desc(parse_fasta(f))
        return dict(seqs)


def parse_greengenes_accessions(f):
    for line in f:
        if line.startswith("#"):
            continue
        line = line.strip()
        yield line.split("\t")
