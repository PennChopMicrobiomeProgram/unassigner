from cStringIO import StringIO

def parse_fasta(f):
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
            seq.write(line)
    yield desc, seq.getvalue()
