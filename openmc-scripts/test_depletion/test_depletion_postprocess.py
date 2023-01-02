import openmc
from openmc.deplete import Results
import sys
import os
from arc_nonproliferation.postprocess import *

if sys.argv[1] is not None:
    base_dir = './' + sys.argv[1]
else:
    raise ValueError("No base directory specified!")

os.chdir(base_dir)

results = Results("depletion_results.h5")
time_steps = results.get_times()
num_steps = len(time_steps)
 
# ==============================================================================
# Fissile Inventory
# ==============================================================================
