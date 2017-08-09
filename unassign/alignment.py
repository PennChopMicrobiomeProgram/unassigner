import shutil
import tempfile

import unassign.parse

class Aligner(object):
    alignment_cls = None
    executable = None

    def __init__(self, species_fp):
        self.species_fp = species_fp
        self.species_max_hits = 1
        self.species_input_fp = None
        self.species_output_fp = None
        self.num_cpus = 1

    @classmethod
    def is_installed(cls):
        exe_path = shutil.which(cls.executable)
        return exe_path is not None

    def search_species(self, seqs):
        """Search species typestrains for match to query sequences."""
        if self.species_input_fp:
            input_fp = self.species_input_fp
            with open(input_fp, "w") as f:
                unassign.parse.write_fasta(f, seqs)
        else:
            infile = tempfile.NamedTemporaryFile(mode="w+t", encoding="utf-8")
            unassign.parse.write_fasta(infile, seqs)
            infile.seek(0)
            input_fp = infile.name

        if self.species_output_fp:
            output_fp = self.species_output_fp
        else:
            outfile = tempfile.NamedTemporaryFile()
            output_fp = outfile.name

        self._call(input_fp, self.species_fp, output_fp)
        return self._load(output_fp)

    @classmethod
    def _load(self, output_fp):
        with open(output_fp) as f:
            hits = [self.alignment_cls(x) for x in self._parse(f)]
        return hits

    @classmethod
    def _parse(cls, f):
        for line in f:
            line = line.strip()
            if line.startswith("#") or (not line):
                continue
            vals = line.split("\t")
            vals = [fn(v) for fn, v in zip(cls.field_types, vals)]
            yield dict(zip(cls.field_names, vals))



class Alignment(object):
    def alignment_length(self):
        raise NotImplemented()

    def num_matches(self):
        raise NotImplemented()

    def unaligned_length(self):
        raise NotImplemented()
