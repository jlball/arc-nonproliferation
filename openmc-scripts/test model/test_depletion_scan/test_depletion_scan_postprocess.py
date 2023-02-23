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

    """ Extract time to 1 SQ for Uranium """
    os.chdir(base_dir + "/Uranium/" + str(mass))

    U_results = Results('depletion_results.h5')
    U_time_to_SQ[i] = extract_time_to_sq('U', U_results)

    os.chdir("../../..")

    """ Extract time to 1 SQ for Thorium """
    os.chdir(base_dir + "/Thorium/" + str(mass))

    Th_results = Results('depletion_results.h5')
    Th_time_to_SQ[i] = extract_time_to_sq('Th', Th_results)

    os.chdir("../../..")

print("Uranium times to 1 SQ:", U_time_to_SQ)
print("Thorium times to 1 SQ:", Th_time_to_SQ)

# ====================================================
# Fission Power
# ====================================================


# ====================================================
# Decay Heat
# ====================================================


# ====================================================
# Plotting
# ====================================================

#Change into dedicated directory for figures or create figures directory
try:
    os.chdir(base_dir + "/figures")
except:
    os.mkdir(base_dir + "/figures")
    os.chdir(base_dir + "/figures")

""" Time to 1 Significant Quantity """
fig, ax = plt.subplots()
ax.spines["top"].set_color("None")
ax.spines["right"].set_color("None")

ax.scatter(masses/1e3, U_time_to_SQ/24, label="$^{238}$U", marker='o')
ax.scatter(masses/1e3, Th_time_to_SQ/24, label="$^{232}$Th", marker='s')
ax.legend()

ax.set_title("Time to 1 Significant Quantity")
ax.set_ylabel("Time (days)")
ax.set_xlabel("Mass of Fertile Material (metric tons)")

fig.savefig("time_to_sq.png")