FROM condaforge/mambaforge:latest

# Setup
WORKDIR /unassigner

COPY unassigner/ /unassigner/unassigner/
COPY README.md requirements.txt setup.py /unassigner/

# Install environment
RUN mamba create --name unassigner -c conda-forge -c bioconda vsearch

ENV PATH="/opt/conda/envs/unassigner/bin/:${PATH}"

# "Activate" the environment
SHELL ["conda", "run", "--no-capture-output", "-n", "unassigner", "/bin/bash", "-c"]

RUN pip install /unassigner/

RUN echo "Python: $(python --version), Conda: $(conda --version), Vsearch: $(vsearch --help)" > installed_packages.txt

# Run
CMD "bash"