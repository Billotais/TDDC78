#!/bin/bash

#SBATCH -J test1
#SBATCH -t 00:10:00
#SBATCH -n 16
#SBATCH -A snic2018-7-5

module load python/2.7.6 
module load buildenv-intel/2015-1

python eval.py
