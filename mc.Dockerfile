FROM python:3.5
RUN apt-get update && apt-get install -yq --no-install-recommends \
    vim \
    postgresql \
    rsync \
    wget \
    docker-io \
    && apt-get clean && \
    rm -rf /var/lib/apt/lists/*
RUN mkdir /mc
WORKDIR /mc/dj

# conda setup
RUN ["/bin/bash", "-c", "wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh"]
RUN chmod 0755 /tmp/miniconda.sh
RUN ["/bin/bash", "-c", "/tmp/miniconda.sh -b -p /conda"]
ENV PATH="/conda/bin:$PATH"
RUN rm /tmp/miniconda.sh
RUN echo 'export PATH=/conda/bin:$PATH' >> /root/.bashrc
RUN conda install --yes conda-build python=3.5
RUN conda install --yes --channel https://conda.anaconda.org/rdkit rdkit
RUN conda install --yes --channel https://conda.anaconda.org/openbabel openbabel

ADD requirements.txt /requirements.txt
RUN pip install -r /requirements.txt
