#!/bin/bash

# Configure the resources required
#SBATCH -p a100                                                # partition (this is the queue your job will be added to)
#SBATCH -N 1
#SBATCH -n 1              	                                # number of tasks (sequential job starts 1 task) (check this if your job unexpectedly uses 2 nodes)
#SBATCH -c 16              	                                # number of cores (sequential job calls a multi-thread program that uses 8 cores)
#SBATCH --time=2-00:00:00
#SBATCH --gres=gpu:1 # generic resource required (here requires 1 GPU)
#SBATCH --mem=100GB                                              # specify memory required per node (here set to 16 GB)
#SBATCH --output=log/02_04_v02.log


#SBATCH --mail-type=ALL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=a1768536@adelaide.edu.au # Where to send mail


# Execute the Python script with srun
date
python -m spanet.train -of options_files/full_hadronic_ttbar/full_training_slurmTest.json --gpus 1
date

