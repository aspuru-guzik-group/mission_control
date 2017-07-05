#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

python -m mc.job_engines.cli \
  run_job_submission \
  $@
