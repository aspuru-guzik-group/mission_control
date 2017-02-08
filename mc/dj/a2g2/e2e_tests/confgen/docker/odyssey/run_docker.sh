docker run \
  -it \
  --rm \
  -h ernie \
  -p 2222:22 \
  -v $(pwd)/docker-entrypoint.sh:/docker-entrypoint.sh \
  -v $(pwd)/slurm.conf:/etc/slurm/slurm.conf \
  -v $(pwd)/slurmdbd.conf:/etc/slurm/slurmdbd.conf \
  -v $(pwd)/supervisord.conf:/etc/supervisord.conf \
  -v ~/anaconda:/anaconda \
  adorsk/slurm
