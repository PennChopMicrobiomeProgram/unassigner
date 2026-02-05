# Unassigner

<!-- Begin badges -->
[![Tests](https://github.com/PennChopMicrobiomeProgram/unassigner/actions/workflows/pr.yml/badge.svg)](https://github.com/PennChopMicrobiomeProgram/unassigner/actions/workflows/pr.yml)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/9ee8fc7bc3e940bb812b35006e95937d)](https://app.codacy.com/gh/PennChopMicrobiomeProgram/unassigner/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![codecov](https://codecov.io/gh/PennChopMicrobiomeProgram/unassigner/graph/badge.svg?token=LAFU84K088)](https://codecov.io/gh/PennChopMicrobiomeProgram/unassigner)
[![PyPI](https://badge.fury.io/py/unassigner.svg)](https://pypi.org/project/unassigner/)
[![Bioconda](https://anaconda.org/bioconda/unassigner/badges/downloads.svg)](https://anaconda.org/bioconda/unassigner/)
[![DockerHub](https://img.shields.io/docker/pulls/chopmicrobiome/unassigner)](https://hub.docker.com/repository/docker/chopmicrobiome/unassigner/)
<!-- End badges -->

Evaluate consistency with named bacterial species for short 16S rRNA
marker gene sequences.

## Summary

The [16S rRNA gene](https://en.wikipedia.org/wiki/16S_ribosomal_RNA)
is found in all bacteria, and its gene sequence is highly
conserved. Amplification and sequencing of bacterial 16S rRNA genes is
a common method used to survey bacterial communities in
[microbiome](https://en.wikipedia.org/wiki/Microbiota)
research. However, high throughput instruments are unable to sequence
the entire gene. Therefore, a short region of the gene is selected for
amplification and sequencing.

The resultant sequences, spanning part of the 16S gene, can be used to
identify the types of bacteria present in a specimen. For example, one
sequence might be assigned to the *Streptococcus* genus based on
sequence similarity. Many programs are available to carry out such
taxonomic assignment.

It is generally thought that the 16S rRNA gene is not suitable for
assignment of bacterial species. We agree, but with a catch: the gene
sequence is suitable for **ruling out** assignment to many bacterial
species. This software is designed to rule out all the species
designations that are inconsistent with a partial 16S rRNA gene
sequence. For those species that are not definitively ruled out, we
assign a probability that the sequence is inconsistent with the
species.

Because the software is geared towards ruling out species rather than
deciding on the best assignment, we call it the *unassigner*. It's a
cheesy joke, but we've decided to roll with it.

The unassigner library provides a command-line program, `unassign`,
that takes a FASTA file of DNA sequences in a 16S gene region, and
gives the probability that the sequence is inconsistent with nearby
bacterial species.

## Installation

Install with conda using:

```bash
conda create --name unassigner -c conda-forge -c bioconda unassigner
```

Or run with Docker using:

```bash
docker run --rm -it chopmicrobiome/unassigner:latest unassign --help
```

### Alternative Installation

Unassigner can be installed using pip:

```bash
pip install unassigner
```

But will require ``vsearch`` to be installed separately. It can also be installed from GitHub:

```bash
git clone https://github.com/PennChopMicrobiomeProgram/unassigner.git
pip install unassigner/
```

## Usage

The `unassign` program requires one argument, a FASTA-formatted file
of short 16S sequences:

```bash
unassign my_sequences.fasta
```

If the program has not been run before, it will automatically download
the bacterial species data it needs, format its reference files,
create an output directory named `my_sequences_unassigned`, and write
a table of results there, along with some auxiliary output files. Note 
that the output directory will be in the same directory as `my_sequences.fasta`.

Please see the output of `unassign --help` for a list of the available
options.

### Trim ragged

The `trimragged` program takes in a query sequence to search and trim and an input fasta file (or it can read from stdin):

```bash
trimragged AGAGTTTGATCCTGGCTCAG --input_file my_sequences.fasta
```

Trimragged is included to extract different regions from the full length 16S rRNA gene. The purpose of this auxiliary software is to account for the full length 16S rRNA sequences where only a part of the primer is present in the sequence. This can be due to low quality at the beginning or at the end of a sequence due to limitations of sequencing platforms. 

The software operates in three steps: 1) Matching the full length of the primer, 2) Matching the partial primer, 3) Aligning reads to other sequences with a known primer location. The sequence of the primer to search and trim is required for the software. Only one primer is accepted at a time, so the user needs to run the software twice with each primer sequence.

Step 1: The software first searches for the full length of the primer sequence. If mismatches are allowed, then the software expands all possibilities of the primer sequence mutations in a list and searches for each. Once a hit is found, the start and end index is stored as a PrimerMatch object.

Step 2: If the min_partial argument is greater than 0, the software then searches for partial matches of the primer in the remaining sequences. The software makes a list of all the possibilities of primers, removing nucleotides from the beginning of the sequence till the minimum length specified by min_partial is reached. Then the software searches for each of the possible primer sequences. Once a hit is found, the start and end index is stored as a Primer Match object.

Step 3: The last part of the software relies on building a database of the sequences with already identified primer sequences from the previous two steps. Then the rest of the reads are aligned against the database of sequences with known primer locations using vsearch. Once a hit is found, and the positions of the primers are estimated by extending the aligned region.

Please see the output of `trimragged --help` for a list of the available
options.

## Contributing

We welcome ideas from our users about how to improve this
software. Please open an issue if you have a question or would like to
suggest a feature.
