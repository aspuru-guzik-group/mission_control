#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

python -m mc.job_module_utils.cli \
  run_job_submission \
  --cfg_file_path="$DIR/cfg/job_submission_runner_cfg.py" \
  $@
