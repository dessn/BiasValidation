#!/bin/bash
#SBATCH --job-name=bias_validation
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8gb
#SBATCH --time=05:00:00
#SBATCH --output=output.log

cd /project2/rkessler/SURVEYS/DES/USERS/parmstrong/dev/bias_validation
python bias_validation.py "$@"


