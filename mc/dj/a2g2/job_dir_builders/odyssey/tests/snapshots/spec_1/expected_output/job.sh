#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=my_partition
#SBATCH --sbatch_1=sbatch_1_value

env_var_1=env_var_1_value
env_var_2=env_var_2_value

module load module_1

body_value
