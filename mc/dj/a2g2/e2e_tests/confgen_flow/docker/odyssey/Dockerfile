FROM centos:7

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
        psmisc \
        rsync

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

COPY slurm.conf /etc/slurm/slurm.conf
COPY slurmdbd.conf /etc/slurm/slurmdbd.conf
COPY supervisord.conf /etc/

RUN yum -y install \
      openssh-server openssh-clients \
      openssl

COPY ./setup_odyssey.sh /setup_odyssey.sh
RUN chmod a+x /setup_odyssey.sh && /setup_odyssey.sh

COPY ./entrypoint.sh /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]
CMD ["/bin/bash"]

# Setup test user.

ARG TEST_USER_USERNAME=test_user
ARG TEST_USER_PASS=test_pass

RUN useradd $TEST_USER_USERNAME -p `openssl passwd -crypt $TEST_USER_PASS`

COPY ./setup_conda.sh /setup_conda.sh
RUN chmod a+rx /setup_conda.sh && su - $TEST_USER_USERNAME -c  /setup_conda.sh

COPY a2g2_conda_env.yml /conda_env.yml
RUN su - $TEST_USER_USERNAME -c "conda env create -f /conda_env.yml -n a2g2_env"
RUN su - $TEST_USER_USERNAME -c "source activate a2g2_env && conda develop /a2g2_env-packages"

RUN yum -y install pam-devel

COPY fake_2factor /fake_2factor
#RUN /fake_2factor/setup_fake_2factor.sh && rm -r /fake_2factor
RUN /fake_2factor/setup_fake_2factor.sh

RUN yum clean all
