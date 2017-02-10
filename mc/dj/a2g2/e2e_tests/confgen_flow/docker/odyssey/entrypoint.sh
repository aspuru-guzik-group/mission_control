#!/bin/bash

echo "Setting slurm host."
sed -i "s/^ControlMachine=.*/ControlMachine=$(hostname)/" /etc/slurm/slurm.conf

echo "Starting supervisor"
/usr/bin/supervisord --configuration /etc/supervisord.conf

exec "$@"
