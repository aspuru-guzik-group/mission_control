#!/bin/bash
#SBATCH --sbatch_1=sbatch_1_value

output_status_file () {
    if [ $? -eq 0 ]; then
        touch ODYSSEY_JOB__COMPLETED
    else
        touch ODYSSEY_JOB__FAILED
    fi
}
trap "output_status_file" EXIT

env_var_1=env_var_1_value
env_var_2=env_var_2_value

module load module_1

body_value
