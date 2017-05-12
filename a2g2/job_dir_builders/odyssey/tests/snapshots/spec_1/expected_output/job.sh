#!/bin/bash
#SBATCH --sbatch_1=sbatch_1_value
#SBATCH --error=ODYSSEY_JOB.stderr
#SBATCH --output=ODYSSEY_JOB.stdout

START_DIR=$PWD
output_status_file () {
    PREV_RETURN_CODE=$?
    pushd $START_DIR > /dev/null
    if [ $PREV_RETURN_CODE -eq 0 ]; then
        touch ODYSSEY_JOB__COMPLETED
    else
        touch ODYSSEY_JOB__FAILED
        echo "tail -n 50 ODYSSEY_JOB.stdout:" >> ODYSSEY_JOB__FAILED
        tail -n 50 ODYSSEY_JOB.stdout >> ODYSSEY_JOB__FAILED
        echo "tail -n 50 ODYSSEY_JOB.stderr:" >> ODYSSEY_JOB__FAILED
        tail -n 50 ODYSSEY_JOB.stderr >> ODYSSEY_JOB__FAILED
        echo "ls -1:" >> ODYSSEY_JOB__FAILED
        ls -1 >> ODYSSEY_JOB__FAILED
    fi
    popd > /dev/null
}
trap "output_status_file" EXIT

env_var_1=env_var_1_value
env_var_2=env_var_2_value

module load module_1

entrypoint_body_value
