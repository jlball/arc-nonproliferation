#/bin/bash

OPENMC_CROSS_SECTIONS="/usr/local/share/xs_data/endfb80_hdf5/cross_sections.xml"
OMP_NUM_THREADS=32

export OPENMC_CROSS_SECTIONS
export OMP_NUM_THREADS

python3 depletion-scan-arc-1-Li6.py Li6_1 1.0
python3 depletion-scan-arc-1-Li6.py Li6_5 5.0
python3 depletion-scan-arc-1-Li6.py Li6_10 10.0
python3 depletion-scan-arc-1-Li6.py Li6_15 15.0
python3 depletion-scan-arc-1-Li6.py Li6_20 20.0
python3 depletion-scan-arc-1-Li6.py Li6_30 30.0
python3 depletion-scan-arc-1-Li6.py Li6_40 40.0
python3 depletion-scan-arc-1-Li6.py Li6_50 50.0