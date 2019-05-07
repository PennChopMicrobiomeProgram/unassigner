import argparse
import logging
import os

from unassigner.algorithm import (
    UnassignAligner, FileAligner, ThresholdAlgorithm,
)
from unassigner.parse import parse_fasta, parse_species_names
from unassigner.prepare_strain_data import download_type_strain_data

def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("query_fasta", type=argparse.FileType("r"),
        help="Query sequences FASTA file")
    p.add_argument("--output_dir",
        help=(
            "Output directory (default: basename of query sequences FASTA "
            "file, plus '_unassigned')"))
    p.add_argument("--type_strain_fasta", default="unassigner_species.fasta",
        help=(
            "Type strain sequences FASTA file (default: %(default)s). "
            "If the default file is not found, sequences are downloaded "
            "and re-formatted automatically."))
    p.add_argument("--num_cpus", type=int,
        help=(
            "Number of CPUs to use during sequence aligment (default: "
            "use all the CPUs)"))
    p.add_argument("--verbose", action="store_true",
        help= "Activate verbose mode.")
    args = p.parse_args(argv)

    if args.verbose is True:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    query_seqs = list(parse_fasta(args.query_fasta, trim_desc=True))

    if args.output_dir is None:
        output_dir = os.path.splitext(args.query_fasta.name)[0] + "_unassigned"
    else:
        output_dir = args.output_dir

    # Download type strain files if needed
    type_strain_fp_is_default = (
        args.type_strain_fasta == p.get_default("type_strain_fasta"))
    type_strain_fp_is_missing = not os.path.exists(args.type_strain_fasta)
    if type_strain_fp_is_default and type_strain_fp_is_missing:
        download_type_strain_data()

    with open(args.type_strain_fasta) as f:
        species_names = dict(parse_species_names(f))

    writer = OutputWriter(output_dir, species_names)

    alignment_query_fp = writer.output_fp("unassigner_query.fasta")
    alignment_output_fp = writer.output_fp("unassigner_query_hits.txt")
    if os.path.exists(alignment_output_fp):
        a = FileAligner(args.type_strain_fasta, alignment_output_fp)
    else:
        a = UnassignAligner(args.type_strain_fasta)
        a.species_input_fp = alignment_query_fp
        a.species_output_fp = alignment_output_fp
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
