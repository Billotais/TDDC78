#!/bin/bash

#SBATCH -J test1
#SBATCH -t 00:55:00
#SBATCH -N 1
#SBATCH -A snic2018-7-5

module load python/2.7.6 
module load buildenv-intel/2015-1

make laplsolv

python eval_lab3.py
