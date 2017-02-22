FROM centos:7

RUN yum makecache fast \
    && yum -y install epel-release \
    && yum -y install \
      bzip2 \
      gcc \
      gcc-c++ \
      make \
      supervisor

RUN yum clean all

ARG CONDA_ROOT=/conda
ARG COMPOSE_ROOT
COPY $COMPOSE_ROOT/setup_conda.sh /setup_conda.sh
ARG MINICONDA_URL
RUN /bin/bash -c "CONDA_ROOT=$CONDA_ROOT MINICONDA_URL=$MINICONDA_URL /setup_conda.sh"

ARG CTX_ROOT
ARG CONDA_ENV_FILE
ARG IMG_CONDA_ENV_FILE=/conda_env.yml
COPY $COMPOSE_ROOT/$CONDA_ENV_FILE $IMG_CONDA_ENV_FILE
RUN /bin/bash -c "source $CONDA_ROOT/bin/activate root && conda env update -f $IMG_CONDA_ENV_FILE --name=root"

COPY $CTX_ROOT/entrypoint.sh /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]

COPY $CTX_ROOT/wait_for_postgres.sh /wait_for_postgres.sh
COPY $CTX_ROOT/start_mc_server.sh /start_mc_server.sh
ENV CONDA_ROOT=$CONDA_ROOT
CMD /bin/bash -c "source \$CONDA_ROOT/bin/activate root && /wait_for_postgres.sh /start_mc_server.sh"
