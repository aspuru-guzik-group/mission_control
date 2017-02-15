FROM ubuntu:yakkety
RUN apt-get update && apt-get install -yq --no-install-recommends \
    bzip2 \
    gcc \
    g++ \
    make \
    curl \
    ca-certificates \
    vim \
    postgresql \
    libpq-dev \
    rsync

# Install docker. We use docker w/in docker to run integration tests.
RUN curl https://apt.dockerproject.org/repo/pool/main/d/docker-engine/docker-engine_1.13.1-0~ubuntu-yakkety_amd64.deb > /tmp/docker-engine.deb
RUN dpkg -i /tmp/docker-engine.deb 2>/dev/null; apt-get -f install -y && dpkg -i /tmp/docker-engine.deb && rm /tmp/docker-engine.deb
RUN curl -L https://github.com/docker/compose/releases/download/1.11.0/docker-compose-`uname -s`-`uname -m` > /usr/local/bin/docker-compose && chmod +x /usr/local/bin/docker-compose

RUN mkdir /mc
WORKDIR /mc/dj

# conda setup
RUN curl https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
  > /tmp/miniconda.sh
RUN chmod 0755 /tmp/miniconda.sh
RUN /tmp/miniconda.sh -b -p /conda
ENV PATH="/conda/bin:$PATH"
RUN rm /tmp/miniconda.sh
RUN echo 'export PATH=/conda/bin:$PATH' >> /root/.bashrc
RUN conda install --yes conda-build python=3.5
RUN conda install --yes --channel https://conda.anaconda.org/rdkit rdkit
RUN conda install --yes --channel https://conda.anaconda.org/openbabel openbabel

ADD requirements.txt /requirements.txt
RUN pip install -r /requirements.txt

RUN apt-get install -y openssh-client

RUN conda install -c conda-forge xorg-libxrender=0.9.10

RUN apt-get clean && rm -rf /var/lib/apt/lists/*
