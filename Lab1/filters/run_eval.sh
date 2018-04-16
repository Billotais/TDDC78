#!/bin/bash

#SBATCH -J test1
#SBATCH -t 01:00:00
#SBATCH -n 128
#SBATCH -A snic2018-7-5

module load python/2.7.12
module load buildenv-intel/2015-1

python eval.py
