import openmc
from openmc.deplete import Results
from arc_nonproliferation.postprocess import *
from arc_nonproliferation.constants import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mpl
import os
from scipy.optimize import curve_fit

masses = [1, 5, 10, 15, 20, 30, 40, 50]
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
        os.chdir(base_dir + "/Uranium/" + str(round(enrichment, 1)).rstrip('0').rstrip('.'))

        U_results = Results('depletion_results.h5')
        U_time_to_SQ[j] = extract_time_to_sq('U', U_results)

        #While we're here, get the number of depletion steps:
        time_steps = U_results.get_times()
        num_steps = len(time_steps)

        os.chdir("../../..")

        """ Extract time to 1 SQ for Thorium """
        os.chdir(base_dir + "/Thorium/" + str(round(enrichment, 1)).rstrip('0').rstrip('.'))

        Th_results = Results('depletion_results.h5')
        Th_time_to_SQ[j] = extract_time_to_sq('Th', Th_results)

        os.chdir("../../..")

    U_times_to_SQ[i] = U_time_to_SQ
    Th_times_to_SQ[i] = Th_time_to_SQ


# ====================================================
# Plotting
# ====================================================


""" Time to 1 Significant Quantity """
text_offset = 5
cm = mpl.cm.Reds
norm = colors.LogNorm(vmin=3, vmax=110)

fig, ax = plt.subplots()
ax.spines["top"].set_color("None")
ax.spines["right"].set_color("None")

# Fit data to 1/x function:
def fit(x, A, B):
    return (A/x) + B

fit_masses = np.linspace(1, masses[-1], num=100)



for i, enrichment in enumerate(enrichments):
    current_color = cm.__call__(norm(enrichment+10))

    ax.scatter(masses, U_times_to_SQ[:, i]/24, label=str(enrichment), marker='o', color=current_color)
    ax.annotate(str(enrichment)+"%", (masses[-1], U_times_to_SQ[-1, i]/24), color=current_color, textcoords='offset points', xytext=(text_offset, 0))

    U_popt, U_pcov = curve_fit(fit, masses, U_times_to_SQ[:, i]/24)
    ax.plot(fit_masses, fit(fit_masses, *U_popt), color=current_color)

#ax.set_xlim(0, enrichments[-1] + 2)
#ax.set_ylim(0, np.max(Th_time_to_SQ/24) + 5)

ax.set_yscale("log")
ax.set_xlim(0, masses[-1] + 2)

ax.set_title("Time to 1 SQ with a Uranium Dopant and Varying Li-6 Enrichment", fontsize=14)
ax.set_ylabel("Time (days)", fontsize=14)
ax.set_xlabel("Mass of Fertile Material (metric tons)", fontsize=14)

fig.set_size_inches(8, 6)

fig.savefig("time_to_sq_U.png", dpi=300)

cm = mpl.cm.Greens

fig, ax = plt.subplots()
ax.spines["top"].set_color("None")
ax.spines["right"].set_color("None")

for i, enrichment in enumerate(enrichments):
    current_color = cm.__call__(norm(enrichment+10))

    ax.scatter(masses, Th_times_to_SQ[:, i]/24, label=str(enrichment), marker='o', color=current_color)
    ax.annotate(str(enrichment)+"%", (masses[-1], Th_times_to_SQ[-1, i]/24), color=current_color, textcoords='offset points', xytext=(text_offset, -text_offset))

    Th_popt, Th_pcov = curve_fit(fit, masses, Th_times_to_SQ[:, i]/24)
    ax.plot(fit_masses, fit(fit_masses, *Th_popt), color=current_color)

ax.set_xlim(0, masses[-1] + 2)
#ax.set_ylim(0, np.max(Th_time_to_SQ/24) + 5)

ax.set_yscale("log")
#ax.set_xscale("log")

ax.set_title("Time to 1 SQ with a Thorium Dopant and Varying Li-6 Enrichment", fontsize=14)
ax.set_ylabel("Time (days)", fontsize=14)
ax.set_xlabel("Mass of Fertile Material (metric tons)", fontsize=14)

fig.set_size_inches(8, 6)

fig.savefig("time_to_sq_Th.png", dpi=300)