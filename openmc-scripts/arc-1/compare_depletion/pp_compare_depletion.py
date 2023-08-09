import openmc
from openmc.deplete import Results
import matplotlib.pyplot as plt
import numpy as np
import sys
from arc_nonproliferation.postprocess import *
from arc_nonproliferation.constants import *


###
# COUPLED
###
os.chdir("COUPLED")

coupled_results = Results("depletion_results.h5")
coupled_masses = get_masses_from_mats("U", coupled_results)

time_steps = coupled_results.get_times()

os.chdir("..")


###
# INDEPENDENT
###
os.chdir("INDEPENDENT")

independent_results = Results("depletion_results.h5")
independent_masses = get_masses_from_mats('U', coupled_results)

os.chdir('..')


###
# PLOTTING
###

fig, ax = plt.subplots()

ax.scatter(time_steps, coupled_masses, label = "Coupled")
ax.scatter(time_steps, independent_masses, label = "Independent")

ax.set_xlabel("Time (days)")
ax.set_ylabel("Mass of Pu-239 (kg)")
ax.set_title("Comparison of Independent and Coupled Operators")

fig.savefig("depletion_comparison.png", dpi=300)