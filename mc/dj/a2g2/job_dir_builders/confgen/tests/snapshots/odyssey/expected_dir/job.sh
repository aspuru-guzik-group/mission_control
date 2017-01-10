#!/bin/bash




module load conda


source /n/aagfs01/software/envs/a2g2_env/bin/activate
python -m a2g2.conformer_generator \
  --input=./conformer_generator.in.json

