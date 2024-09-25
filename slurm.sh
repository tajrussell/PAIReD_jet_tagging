#!/bin/bash

#SBATCH -p gpu --gres=gpu:2
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=192g
#SBATCH --time=48:00:00
#SBATCH --exclude=gpu2509
#SBATCH -J  Classi
#SBATCH -o  /users/trussel1/PAIReD_jet_tagging/log/run-Classi-%j.out  # File to which STDOUT will be written
#SBATCH -e  /users/trussel1/PAIReD_jet_tagging/log/run-Classi-%j.out  # File to which STDERR will be written
#SBATCH --mail-type FAIL                 # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user trevor_russell@brown.edu  # Email to which notifications will be sent

INDIR=/users/trussel1/PAIReD_jet_tagging/
export PATH="/users/trussel1/miniconda3/envs/weaver/bin/:$PATH"
cd $INDIR
. train.sh
