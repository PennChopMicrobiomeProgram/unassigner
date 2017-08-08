import shutil


class Aligner(object):
    alignment_cls = None
    executable = None

    def __init__(self, species_fp):
        self.species_fp = species_fp

    @classmethod
    def is_installed(cls):
        exe_path = shutil.which(cls.executable)
        return exe_path is not None

    def search_species(self, seqs):
        raise NotImplemented()


class Alignment(object):
    def alignment_length(self):
        raise NotImplemented()

    def num_matches(self):
        raise NotImplemented()

    def unaligned_length(self):
        raise NotImplemented()
