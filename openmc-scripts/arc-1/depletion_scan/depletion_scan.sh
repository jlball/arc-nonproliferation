#!/bin/bash
#SBATCH -n 128
#SBATCH -N 8
#SBATCH -t 0-01:00
#SBATCH -p sched_mit_hill
#SBATCH --mem-per-cpu=4000
#SBATCH -o output_%j.txt
#SBATCH -e error_%j.txt
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=jlball@mit.edu
module load gcc/11.2.0
module load cmake/3.17.3
module load openmpi/4.1.2
module load hdf5/1.10.5-parallel
module load python/3.9.4
pip install --user "/home/jlball/openmc"
export OPENMC_CROSS_SECTIONS="/home/jlball/xs_data/endfb-viii.0-hdf5/cross_sections.xml"
python depletion-scan-arc-1.py engaging_run