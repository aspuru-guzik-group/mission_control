FROM centos:7

RUN yum makecache fast \
    && yum -y install epel-release \
    && yum -y install \
      bzip2 \
      gcc \
      gcc-c++ \
      make \
      postgresql \
      postgresql-devel \
      supervisor

RUN yum clean all

ARG CONDA_ROOT=/root/conda
COPY ./setup_conda.sh /setup_conda.sh
RUN chmod a+x /setup_conda.sh && CONDA_ROOT=$CONDA_ROOT /setup_conda.sh

COPY ./entrypoint.sh /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]

ARG CONDA_ENV_FILENAME=mc_conda_env.yml
COPY ./$CONDA_ENV_FILENAME /$CONDA_ENV_FILENAME
RUN $CONDA_ROOT/bin/conda env update -f /$CONDA_ENV_FILENAME

COPY ./start_mc_server.sh /start_mc_server.sh
CMD /start_mc_server.sh
