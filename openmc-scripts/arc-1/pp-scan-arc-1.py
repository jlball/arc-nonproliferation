import openmc
import matplotlib.pyplot as plt
import numpy as np
import sys
from arc_nonproliferation.postprocess import *
from arc_nonproliferation.constants import *

if sys.argv[1] is not None:
    base_dir = './' + sys.argv[1]
else:
    raise ValueError("No base directory specified!")

fusion_power = 500 #MW
total_neutron_rate = fusion_power * neutrons_per_MJ

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

    """ Extract time to 1 SQ for Uranium """
    os.chdir(base_dir + "/Uranium/" + str(mass))
    sp = openmc.StatePoint('statepoint.10.h5')
    U_tally = sp.get_tally(name='Breeding Tally')
    breeding_rate = U_tally.get_values(scores=['absorption'], nuclides=['U238'])
    print(breeding_rate)
    
    #U_time_to_SQ[i]

    os.chdir("../../..")

    # """ Extract time to 1 SQ for Thorium """
    # os.chdir(base_dir + "/Thorium/" + str(mass))

    # Th_results = Results('depletion_results.h5')
    # Th_time_to_SQ[i] = extract_time_to_sq('Th', Th_results)

    # os.chdir("../../..")
