# Unassigner

[![Tests](https://github.com/Ulthran/unassigner/actions/workflows/tests.yml/badge.svg)](https://github.com/Ulthran/unassigner/actions/workflows/tests.yml)
[![Super-Linter](https://github.com/Ulthran/unassigner/actions/workflows/linter.yml/badge.svg)](https://github.com/Ulthran/unassigner/actions/workflows/linter.yml)
[![codecov](https://codecov.io/gh/Ulthran/unassigner/branch/master/graph/badge.svg?token=76YWJFWGON)](https://codecov.io/gh/Ulthran/unassigner)
[![PyPi](https://github.com/Ulthran/unassigner/actions/workflows/python-publish.yml/badge.svg)](https://github.com/Ulthran/unassigner/actions/workflows/python-publish.yml)

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

The Python library and command-line program can be installed using
[pip](https://pypi.org/project/pip/).

```bash
pip install unassigner
```

Besides the python libraries listed in the setup file, this program
requires `vsearch` to be installed.  This program is used to search
for the closest matching bacterial species and return pairwise
sequence alignments.  It's available through
[conda](https://anaconda.org/bioconda/vsearch), and this is our
recommended method for installation.

```bash
conda create --name unassigner
conda activate unassigner
conda install -c bioconda vsearch
pip install unassigner
```

### Alternative Installation

If `pip install unassigner` isn't working or if you want to use a development 
version, you can also install via git.

```bash
conda create --name unassigner
conda activate unassigner
conda install -c bioconda vsearch
git clone https://github.com/kylebittinger/unassigner.git
cd unassigner
pip install -r requirements.txt
pip install .
```

If you don't want to use conda, see the 
[vsearch repo](https://github.com/torognes/vsearch) for alternative install 
methods.

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

## Contributing

We welcome ideas from our users about how to improve this
software. Please open an issue if you have a question or would like to
suggest a feature.
