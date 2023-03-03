import openmc
import matplotlib.pyplot as plt
import numpy as np
import sys
from arc_nonproliferation.postprocess import *
from arc_nonproliferation.constants import *
from scipy.optimize import curve_fit
from scipy.stats import linregress

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
enrichments = np.loadtxt(base_dir + '/enrichments.txt')

U_time_to_SQ = np.empty(len(enrichments))
Th_time_to_SQ = np.empty(len(enrichments))

""" Iterate through each mass simulated and compute time to SQ"""
for i, enrichment in enumerate(enrichments):

    """ Extract time to 1 SQ for Uranium """
    os.chdir(base_dir + "/Uranium/" + str(enrichment))
    sp = openmc.StatePoint('statepoint.10.h5')
    U_tally = sp.get_tally(name='Breeding Tally')

    breeding_rate = U_tally.get_values(scores=['absorption'], nuclides=['U238']) * total_neutron_rate * Pu239_mass_in_kg
    U_time_to_SQ[i] = sig_quantity / breeding_rate # seconds
    os.chdir("../../..")

    """ Extract time to 1 SQ for Thorium """
    os.chdir(base_dir + "/Thorium/" + str(enrichment))
    sp = openmc.StatePoint('statepoint.10.h5')
    Th_tally = sp.get_tally(name='Breeding Tally')

    breeding_rate = Th_tally.get_values(scores=['absorption'], nuclides=['Th232']) * total_neutron_rate * U233_mass_in_kg
    Th_time_to_SQ[i] = sig_quantity / breeding_rate # seconds
    os.chdir("../../..")

#Convert to units of hours
U_time_to_SQ = U_time_to_SQ/3600
Th_time_to_SQ = Th_time_to_SQ/3600

# ====================================================
# TBR
# ====================================================

U_tbr = np.empty(len(enrichments))
Th_tbr = np.empty(len(enrichments))

for i, enrichment in enumerate(enrichments):

    """ Extract time to 1 SQ for Uranium """
    os.chdir(base_dir + "/Uranium/" + str(enrichment))
    sp = openmc.StatePoint('statepoint.10.h5')
    U_tally = sp.get_tally(name='FLiBe Tally')

    U_tbr[i] = U_tally.get_values(scores=['(n,Xt)'])
    os.chdir("../../..")

    """ Extract time to 1 SQ for Thorium """
    os.chdir(base_dir + "/Thorium/" + str(enrichment))
    sp = openmc.StatePoint('statepoint.10.h5')
    Th_tally = sp.get_tally(name='FLiBe Tally')

    Th_tbr[i] = Th_tally.get_values(scores=['(n,Xt)'])
    os.chdir("../../..")

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

ax.scatter(enrichments, U_time_to_SQ/24, label="$^{238}$U", marker='o')
ax.scatter(enrichments, Th_time_to_SQ/24, label="$^{232}$Th", marker='s')

np.save("U_time_to_SQ_normal", U_time_to_SQ)
np.save("Th_time_to_SQ_normal", Th_time_to_SQ)

ax.legend()

ax.set_xlim(0, enrichments[-1] + 2)
ax.set_ylim(0, np.max(Th_time_to_SQ/24) + 5)

ax.set_title("Time to Breed a Significant Quantity of Fissile Material", fontsize=14)
ax.set_ylabel("Time (days)", fontsize=14)
ax.set_xlabel("Li6 Enrichment", fontsize=14)

fig.savefig("time_to_sq.png")

# TBR
fig, ax = plt.subplots()
ax.spines["top"].set_color("None")
ax.spines["right"].set_color("None")

ax.scatter(enrichments, U_tbr, label="$^{238}$U", marker='o')
ax.scatter(enrichments, Th_tbr, label="$^{232}$Th", marker='s')

ax.legend()

ax.set_title("Tritium Breeding Ratio", fontsize=14)
ax.set_ylabel("TBR", fontsize=14)
ax.set_xlabel("Li6 Enrichment", fontsize=14)

ax.set_xlim(0, enrichments[-1] + 2)
ax.set_ylim(0.6, 1.2)

fig.savefig("tbr.png")