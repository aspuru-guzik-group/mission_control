FROM ubuntu:yakkety
RUN apt-get update && apt-get install -yq --no-install-recommends \
    bzip2 \
    gcc \
    g++ \
    make \
    curl \
    ca-certificates \
    vim \
    rsync

# Install docker. We use docker w/in docker to run integration tests.
RUN curl https://apt.dockerproject.org/repo/pool/main/d/docker-engine/docker-engine_1.13.1-0~ubuntu-yakkety_amd64.deb > /tmp/docker-engine.deb
RUN dpkg -i /tmp/docker-engine.deb 2>/dev/null; apt-get -f install -y && dpkg -i /tmp/docker-engine.deb && rm /tmp/docker-engine.deb
RUN curl -L https://github.com/docker/compose/releases/download/1.11.0/docker-compose-`uname -s`-`uname -m` > /usr/local/bin/docker-compose && chmod +x /usr/local/bin/docker-compose

RUN mkdir /mc
WORKDIR /mc/dj

# conda setup
ARG MINICONDA_URL
RUN curl $MINICONDA_URL > /tmp/miniconda.sh
RUN chmod 0755 /tmp/miniconda.sh
ARG CONDA_ROOT=/conda/
RUN /tmp/miniconda.sh -b -p $CONDA_ROOT
RUN rm /tmp/miniconda.sh
RUN echo "export PATH=$CONDA_ROOT/bin:\$PATH" >> /root/.bashrc
RUN $CONDA_ROOT/bin/conda install --yes conda-build
ARG CONDA_ENV_FILE=a2g2_mc_conda_env.yml
COPY $CONDA_ENV_FILE /$CONDA_ENV_FILE
RUN $CONDA_ROOT/bin/conda env update -f /$CONDA_ENV_FILE --name=root
#ARG PIP_REQUIREMENTS_FILE=a2g2_mc_pip_requirements.txt
#COPY $PIP_REQUIREMENTS_FILE /$PIP_REQUIREMENTS_FILE
#RUN /bin/bash -c "source $CONDA_ROOT/bin/activate root && pip install -r /$PIP_REQUIREMENTS_FILE"

RUN apt-get install -y openssh-client

RUN apt-get clean && rm -rf /var/lib/apt/lists/*
