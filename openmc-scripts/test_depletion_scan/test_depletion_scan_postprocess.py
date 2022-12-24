import openmc
from openmc.deplete import Results
import matplotlib.pyplot as plt
import numpy as np
import sys
from arc_nonproliferation.postprocess import *

if sys.argv[1] is not None:
    base_dir = './' + sys.argv[1]
else:
    raise ValueError("No base directory specified!")

# ====================================================
# Time to a Significant Quantity
# ====================================================

""" Load masses and initialisze final output arrays """
masses = np.loadtxt(base_dir + '/masses.txt')
print(masses)
U_time_to_SQ = np.empty(len(masses))
Th_time_to_SQ = np.empty(len(masses))

""" Iterate through each mass simulated and compute time to SQ"""
for i, mass in enumerate(masses):
    """ Load depletion results file """
    os.chdir(base_dir + "/Uranium/" + str(mass))

    U_results = Results('depletion_results.h5')
    U_time_to_SQ[i] = extract_time_to_sq('U', U_results)

    os.chdir("../../..")