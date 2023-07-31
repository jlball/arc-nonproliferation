#!/bin/bash
#SBATCH -n 128
#SBATCH -N 4
#SBATCH -t 0-01:00
#SBATCH -p sched_mit_nse
#SBATCH --mem-per-cpu=4000
#SBATCH -o output_%j.txt
#SBATCH -e error_%j.txt
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=jlball@mit.edu
module use /home/software/psfc/modulefiles
module load psfc/openmc
module load python/3.9.4
export OPENMC_CROSS_SECTIONS="/home/jlball/xs_data/endfb-viii.0-hdf5/cross_sections.xml"
cd /home/jlball/arc-nonproliferation
pip3 install --user -e .

cd openmc-scripts/arc-1/depletion
mkdir test_run

export OMP_NUM_THREADS=32
mpirun -n 4 python3 depletion-arc-1-HPC.py 50000 U