#!/bin/bash

set -o errexit

echo "Setting slurm host."
sed -i "s/^ControlMachine=.*/ControlMachine=$(hostname)/" /etc/slurm/slurm.conf

echo "Starting supervisor"
/usr/bin/supervisord --configuration /etc/supervisord.conf

echo "Setting symlink for faked a2g2_env"
ODYSSEY_A2G2_CONDA_ENV_PATH="${ODYSSEY_A2G2_CONDA_ENV_PATH:-/path/to/odyssey_a2g2_conda_env}"
mkdir -p $(dirname $ODYSSEY_A2G2_CONDA_ENV_PATH)
ln -n -s -f ${A2G2_CONDA_ENV_PATH-/path/to/a2g2_conda_env} ${ODYSSEY_A2G2_CONDA_ENV_PATH}

exec "$@"
