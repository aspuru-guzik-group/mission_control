FROM ubuntu:latest

# Install system dependencies.
RUN apt-get update
RUN apt-get install -yq --no-install-recommends openbabel
RUN apt-get install -yq --no-install-recommends wget
RUN apt-get install -yq --no-install-recommends bzip2

# Install conda.
RUN mkdir -p /opt
RUN wget "https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh" \
  -O /tmp/miniconda.sh \
  --no-check-certificate


RUN chmod a+x /tmp/miniconda.sh && /tmp/miniconda.sh -b -p /opt/conda
ENV PATH="/opt/conda/bin:$PATH"

RUN conda install --yes --channel https://conda.anaconda.org/rdkit rdkit python=3

RUN rm /tmp/miniconda.sh
RUN apt-get clean && rm -rf /var/lib/apt/lists/*
