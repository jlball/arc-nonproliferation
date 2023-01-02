import openmc
from openmc.deplete import Results
import sys
import os
from arc_nonproliferation.postprocess import *
import matplotlib.pyplot as plt

if sys.argv[1] is not None:
    base_dir = './' + sys.argv[1]
else:
    raise ValueError("No base directory specified!")

os.chdir(base_dir)

results = Results("depletion_results.h5")
time_steps = results.get_times('h')
num_steps = len(time_steps)
 
# ==============================================================================
# Fissile Inventory
# ==============================================================================
U_fissile_masses = get_masses_from_mats("U", results)


# ==============================================================================
# Plotting
# ==============================================================================

#Change into dedicated directory for figures or create figures directory
try:
    os.chdir("figures")
except:
    os.mkdir("figures")
    os.chdir("figures")

""" Fissile Mass """
fig, ax = plt.subplots()
ax.scatter(time_steps, U_fissile_masses, label="Pu-239")
ax.legend()

ax.set_title("Mass of Pu-239 vs. Time")
ax.set_xlabel("Time (hours)")
ax.set_ylabel("Mass of Fissile Material (kg)")

fig.savefig("fissile_mass.png")
