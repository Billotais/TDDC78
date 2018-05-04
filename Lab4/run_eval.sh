#!/bin/bash

#SBATCH -J test1
#SBATCH -t 01:00:00
#SBATCH -n 16
#SBATCH -A snic2018-7-5

module load python/2.7.6 
module load buildenv-intel/2015-1

make parallel


python eval_lab4.py
