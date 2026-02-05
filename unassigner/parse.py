import logging
from io import StringIO


def parse_species_names(f):
    for desc, _ in parse_fasta(f):
        vals = desc.split("\t", maxsplit=1)
        accession = vals[0]
        if len(vals) == 2:
            species_name = vals[1]
        else:
            species_name = accession
        yield accession, species_name


def parse_fasta(f, trim_desc=False):
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
    try:
        desc = next(f).strip()[1:]
        if trim_desc:
            desc = desc.split()[0]
    except StopIteration:
        return
    seq = StringIO()
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            yield desc, seq.getvalue()
            desc = line[1:]
            if trim_desc:
                desc = desc.split()[0]
            seq = StringIO()
        else:
            seq.write(line.replace(" ", "").replace("U", "T"))
    yield desc, seq.getvalue()


def parse_desc(desc):
    arr = desc.split("|")
    if len(arr) >= 4:
        accession = arr[2]
        species_name = arr[3]
        return accession, species_name

    desc = desc.lstrip(">").strip()
    toks = desc.split()
    if len(toks) >= 2:
        accession = toks[0]
        species_name = " ".join(toks[1:])
        return accession, species_name

    logging.error(f"Couldn't find accession and/or organism identifier in {desc}")
    logging.error("Skipping this sequence...")
    return None, None


def write_fasta(f, seqs):
    for desc, seq in seqs:
        f.write(">{0}\n{1}\n".format(desc, seq))


def load_fasta(filepath, trim_desc=True):
    """Load all sequences from a FASTA file

    Parameters
    ----------
    fasta_fp : Input filepath, FASTA format.

    Returns
    -------
    A dictionary mapping sequence identifiers to sequences.
    """
    with open(filepath) as f:
        seqs = parse_fasta(f, trim_desc=trim_desc)
        return dict(seqs)


def parse_greengenes_accessions(f):
    for line in f:
        if line.startswith("#"):
            continue
        line = line.strip()
        yield line.split("\t")


def cast_num_or_na(val, cast_func):
    if val == "NA":
        return None
    return cast_func(val)


def parse_results(f):
    float_fields = ["probability_incompatible"]
    int_fields = ["region_mismatches", "region_positions"]
    header_line = next(f)
    header_line = header_line.rstrip()
    fields = header_line.split("\t")
    for line in f:
        line = line.rstrip()
        vals = line.split("\t")
        res = dict(zip(fields, vals))
        for field, val in res.items():
            if field in float_fields:
                res[field] = cast_num_or_na(val, float)
            elif field in int_fields:
                res[field] = cast_num_or_na(val, int)
        yield res
