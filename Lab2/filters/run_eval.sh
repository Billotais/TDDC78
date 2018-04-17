#!/bin/bash

#SBATCH -J test1
#SBATCH -t 00:15:00
#SBATCH -n 64
#SBATCH -A snic2018-7-5

module load python/2.7.6 
module load buildenv-intel/2015-1

make blur
make thres

python eval_blur.py
python eval_thres.py
