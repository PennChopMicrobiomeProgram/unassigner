import argparse
import os

from unassign.algorithm import (
    UnassignAligner, FileAligner, ThresholdAlgorithm,
)
from unassign.parse import parse_fasta, parse_species_names

def main(argv=None):
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("query_fasta", type=argparse.FileType("r"),
        help="Query sequences filepath (FASTA format)")
    p.add_argument("output_dir", help="Output directory")
    p.add_argument("--type_strain_fasta", default="species.fasta",
        help="Type strain sequences filepath (FASTA format)")
    p.add_argument("--alignment_file", type=argparse.FileType("r"),
        help="Use pre-computed alignment file")
    p.add_argument("--verbose", action="store_true",
        help= "Activate verbose mode.")
    args = p.parse_args(argv)

    if args.verbose is True:
        logging.basicConfig(level=logging.INFO)

    query_seqs = list(parse_fasta(args.query_fasta, trim_desc=True))

    with open(args.type_strain_fasta) as f:
        species_names = dict(parse_species_names(f))
    writer = OutputWriter(args.output_dir, species_names)

    if args.alignment_file:
        a = FileAligner(args.type_strain_fasta, args.alignment_file)
    else:
        a = UnassignAligner(args.type_strain_fasta)
        a.species_input_fp = writer.output_fp("unassigner_query.fasta")
        a.species_output_fp = writer.output_fp("unassigner_query_blastn.txt")

    algorithm = ThresholdAlgorithm(a)
    for query_id, query_results in algorithm.unassign(query_seqs):
        writer.write_results(query_id, query_results)


class OutputWriter:
    standard_keys = ["typestrain_id", "probability_incompatible"]

    def __init__(self, output_dir, species_names):
        self.species_names = species_names

        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        self.output_dir = output_dir

        standard_fp = self.output_fp("unassigner_output.tsv")
        self.standard_file = open(standard_fp, "w")
        standard_header = ["query_id", "species"] + self.standard_keys
        self._write_tsv_line(self.standard_file, standard_header)

        algorithm_fp = self.output_fp("algorithm_output.tsv")
        self.algorithm_file = open(algorithm_fp, "w")
        self.algorithm_keys = None

    def output_fp(self, filename):
        return os.path.join(self.output_dir, filename)

    @staticmethod
    def _write_tsv_line(fileobj, vals):
        vals = [str(val) for val in vals]
        line = "\t".join(vals)
        fileobj.write(line)
        fileobj.write("\n")

    def write_results(self, query_id, results):
        for res in results:
            # Set custom algorithm keys based on the first result
            if self.algorithm_keys is None:
                self.algorithm_keys = [
                    k for k in res if k not in self.standard_keys]
                algorithm_header = ["query_id"] + self.algorithm_keys
                self._write_tsv_line(self.algorithm_file, algorithm_header)
            species_name = self.species_names.get(res["typestrain_id"])

            standard_vals = [query_id, species_name] + \
                [res[k] for k in self.standard_keys]
            self._write_tsv_line(self.standard_file, standard_vals)
            algorithm_vals = [query_id] + [res[k] for k in self.algorithm_keys]
            self._write_tsv_line(self.algorithm_file, algorithm_vals)
