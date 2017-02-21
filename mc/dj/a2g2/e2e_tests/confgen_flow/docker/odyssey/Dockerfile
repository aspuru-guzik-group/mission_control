FROM centos:7

ARG COMPOSE_ROOT=.
ARG CTX_ROOT=.

ENV SLURM_VERSION 16.05.9
ENV SLURM_DOWNLOAD_MD5 c0e7fe9d9beaf5ca4599e554f1029fc3
ENV SLURM_DOWNLOAD_URL https://www.schedmd.com/downloads/latest/slurm-"$SLURM_VERSION".tar.bz2

RUN yum makecache fast \
    && yum -y install epel-release \
    && yum -y install \
        wget \
        bzip2 \
        gcc \
        gcc-c++\
        make \
        munge \
        munge-devel \
        supervisor \
        mariadb-server \
        mariadb-devel \
        psmisc

RUN groupadd -r slurm && useradd -r -g slurm slurm

RUN set -x \
    && wget -O slurm.tar.bz2 "$SLURM_DOWNLOAD_URL" \
    && echo "$SLURM_DOWNLOAD_MD5" slurm.tar.bz2 | md5sum -c - \
    && mkdir /usr/local/src/slurm \
    && tar jxf slurm.tar.bz2 -C /usr/local/src/slurm --strip-components=1 \
    && rm slurm.tar.bz2 \
    && cd /usr/local/src/slurm \
    && ./configure --enable-debug --enable-front-end --prefix=/usr \
       --sysconfdir=/etc/slurm --with-mysql_config=/usr/bin \
       --libdir=/usr/lib64 \
    && make install \
    && install -D -m644 etc/cgroup.conf.example /etc/slurm/cgroup.conf.example \
    && install -D -m644 etc/slurm.conf.example /etc/slurm/slurm.conf.example \
    && install -D -m644 etc/slurm.epilog.clean /etc/slurm/slurm.epilog.clean \
    && install -D -m644 etc/slurmdbd.conf.example /etc/slurm/slurmdbd.conf.example \
    && install -D -m644 contribs/slurm_completion_help/slurm_completion.sh /etc/profile.d/slurm_completion.sh \
    && cd \
    && rm -rf /usr/local/src/slurm \
    && mkdir /etc/sysconfig/slurm \
        /var/spool/slurmd \
        /var/run/slurmd \
        /var/lib/slurmd \
        /var/log/slurm \
    && /sbin/create-munge-key

COPY $CTX_ROOT/slurm.conf /etc/slurm/slurm.conf
COPY $CTX_ROOT/slurmdbd.conf /etc/slurm/slurmdbd.conf
COPY $CTX_ROOT/supervisord.conf /etc/

RUN yum -y install \
      openssh-server openssh-clients \
      openssl \
      rsync

COPY $CTX_ROOT/setup_odyssey.sh /setup_odyssey.sh
RUN chmod a+x /setup_odyssey.sh && /setup_odyssey.sh

# Fake 2factor auth, via pam/sshd.
RUN yum -y install pam-devel
COPY $CTX_ROOT/fake_2factor /fake_2factor
RUN /fake_2factor/setup_fake_2factor.sh && rm -r /fake_2factor

# Setup test user.
ARG TEST_USER_USERNAME=test_user
ARG TEST_USER_PASS=test_pass

RUN useradd $TEST_USER_USERNAME -p `openssl passwd -crypt $TEST_USER_PASS`

ARG CONDA_ROOT=/conda
COPY $COMPOSE_ROOT/setup_conda.sh /setup_conda.sh
ARG MINICONDA_URL
RUN /bin/bash -c "CONDA_ROOT=$CONDA_ROOT MINICONDA_URL=$MINICONDA_URL /setup_conda.sh"

ARG IMG_CONDA_ENV_FILE=/conda_env.yml
ARG CONDA_ENV_FILE
COPY $COMPOSE_ROOT/$CONDA_ENV_FILE $IMG_CONDA_ENV_FILE
ARG A2G2_CONDA_ENV_PATH="/a2g2_env"
RUN mkdir -p $A2G2_CONDA_ENV_PATH && chown -R test_user $A2G2_CONDA_ENV_PATH
RUN su - $TEST_USER_USERNAME -c "source $CONDA_ROOT/bin/activate && conda env create -f $IMG_CONDA_ENV_FILE --prefix $A2G2_CONDA_ENV_PATH"
RUN su - $TEST_USER_USERNAME -c "source $CONDA_ROOT/bin/activate /$A2G2_CONDA_ENV_PATH && conda develop /a2g2_env-packages"
ENV A2G2_CONDA_ENV_PATH=$A2G2_CONDA_ENV_PATH

RUN yum clean all
