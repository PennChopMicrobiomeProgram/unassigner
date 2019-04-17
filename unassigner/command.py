import argparse
import os

from unassign.algorithm import (
    UnassignAligner, FileAligner, ThresholdAlgorithm,
)
from unassign.parse import parse_fasta, parse_species_names

def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("query_fasta", type=argparse.FileType("r"),
        help="Query sequences FASTA file")
    p.add_argument("--output_dir",
        help="Output directory (default: derived from QUERY_FASTA)")
    p.add_argument("--type_strain_fasta", default="species.fasta",
        help="Type strain sequences FASTA file (default: %(default)s)")
    p.add_argument("--alignment_file", type=argparse.FileType("r"),
        help="Use pre-computed alignment file (default: run new alignments)")
    p.add_argument("--num_cpus", type=int,
        help="Number of CPUs to use in aligment (default: use all the CPUs)")
    p.add_argument("--verbose", action="store_true",
        help= "Activate verbose mode.")
    args = p.parse_args(argv)

    if args.verbose is True:
        logging.basicConfig(level=logging.INFO)

    query_seqs = list(parse_fasta(args.query_fasta, trim_desc=True))

    if args.output_dir is None:
        output_dir = os.path.splitext(args.query_fasta.name)[0] + "_unassigned"
    else:
        output_dir = args.output_dir

    with open(args.type_strain_fasta) as f:
        species_names = dict(parse_species_names(f))

    writer = OutputWriter(output_dir, species_names)

    if args.alignment_file:
        a = FileAligner(args.type_strain_fasta, args.alignment_file)
    else:
        a = UnassignAligner(args.type_strain_fasta)
        a.species_input_fp = writer.output_fp("unassigner_query.fasta")
        a.species_output_fp = writer.output_fp("unassigner_query_hits.txt")
        a.num_cpus = args.num_cpus

    algorithm = ThresholdAlgorithm(a)
    for query_id, query_results in algorithm.unassign(query_seqs):
        writer.write_results(query_id, query_results)


class OutputWriter:
    standard_keys = [
        "typestrain_id", "region_mismatches", "region_positions",
        "probability_incompatible"]

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
        for result in results:
            species_name = self.species_names.get(result["typestrain_id"], "NA")

            # Set custom algorithm keys based on the first result
            if self.algorithm_keys is None:
                self.algorithm_keys = list(result.keys())
                algorithm_header = ["query_id", "species_name"] + self.algorithm_keys
                self._write_tsv_line(self.algorithm_file, algorithm_header)

            standard_vals = [query_id, species_name] + \
                [result[k] for k in self.standard_keys]
            self._write_tsv_line(self.standard_file, standard_vals)

            algorithm_vals = [query_id, species_name] + \
                [result[k] for k in self.algorithm_keys]
            self._write_tsv_line(self.algorithm_file, algorithm_vals)
