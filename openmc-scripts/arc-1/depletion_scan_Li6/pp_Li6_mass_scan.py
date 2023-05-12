import openmc
from openmc.deplete import Results
from arc_nonproliferation.postprocess import *
from arc_nonproliferation.constants import *
import numpy as np
import matplotlib.pyplot as plt
import os

masses = [1, 5, 10, 15, 20, 30]
enrichments = [0, 1, 2, 4, 7.5, 12, 25, 50, 100]
prefix = "Li6_"

directories = []
for i, mass in enumerate(masses):
    directories.append(prefix + str(mass))

U_times_to_SQ = np.empty((len(masses), len(enrichments)))
Th_times_to_SQ = np.empty((len(masses), len(enrichments)))

for i, base_dir in enumerate(directories):

    """ Load masses and initialisze final output arrays """
    enrichments = np.loadtxt(base_dir + '/enrichments.txt')
    U_time_to_SQ = np.empty(len(enrichments))
    Th_time_to_SQ = np.empty(len(enrichments))

    """ Iterate through each mass simulated and compute time to SQ"""
    for j, enrichment in enumerate(enrichments):

        """ Extract time to 1 SQ for Uranium """
        os.chdir(base_dir + "/Uranium/" + str(enrichment))

        U_results = Results('depletion_results.h5')
        U_time_to_SQ[j] = extract_time_to_sq('U', U_results)

        #While we're here, get the number of depletion steps:
        time_steps = U_results.get_times()
        num_steps = len(time_steps)

        os.chdir("../../..")

        """ Extract time to 1 SQ for Thorium """
        os.chdir(base_dir + "/Thorium/" + str(enrichment))

        Th_results = Results('depletion_results.h5')
        Th_time_to_SQ[j] = extract_time_to_sq('Th', Th_results)

        os.chdir("../../..")

    U_times_to_SQ[i] = U_time_to_SQ
    Th_times_to_SQ[i] = Th_time_to_SQ

print(U_times_to_SQ)
