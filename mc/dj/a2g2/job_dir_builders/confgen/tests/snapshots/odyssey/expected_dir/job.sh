#!/bin/bash

module load conda

source activate /n/aagfs1/a2g2
python -m a2g2.conformer_generator --input=./conformer_generator.in.json
