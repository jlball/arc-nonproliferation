import openmc
from openmc.deplete import Results
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
masses = np.loadtxt(base_dir + '/masses.txt')

U_time_to_SQ = np.empty(len(masses))
Th_time_to_SQ = np.empty(len(masses))

""" Iterate through each mass simulated and compute time to SQ"""
for i, mass in enumerate(masses):

    """ Extract time to 1 SQ for Uranium """
    os.chdir(base_dir + "/Uranium/" + str(mass))

    U_results = Results('depletion_results.h5')
    U_time_to_SQ[i] = extract_time_to_sq('U', U_results)

    #While we're here, get the number of depletion steps:
    time_steps = U_results.get_times()
    num_steps = len(time_steps)

    os.chdir("../../..")

    """ Extract time to 1 SQ for Thorium """
    os.chdir(base_dir + "/Thorium/" + str(mass))

    Th_results = Results('depletion_results.h5')
    Th_time_to_SQ[i] = extract_time_to_sq('Th', Th_results)

    os.chdir("../../..")

""" Iterate through each mass simulated and get fissile mass at each time step"""
U_fissile_masses = np.empty((len(masses), len(time_steps)))
Th_fissile_masses = np.empty((len(masses), len(time_steps)))
for i, mass in enumerate(masses):

    os.chdir(base_dir + "/Uranium/" + str(mass))

    U_results = Results('depletion_results.h5')
    U_fissile_masses[i] = get_masses_from_mats('U', U_results)

    os.chdir("../../..")

    os.chdir(base_dir + "/Thorium/" + str(mass))

    Th_results = Results('depletion_results.h5')
    Th_fissile_masses[i] = get_masses_from_mats('Th', Th_results)

    os.chdir("../../..")

# ====================================================
# Plotting
# ====================================================

masses = masses / 1e3 #convert from kg to metric tons

#Change into dedicated directory for figures or create figures directory
try:
    os.chdir(base_dir + "/figures")
except:
    os.mkdir(base_dir + "/figures")
    os.chdir(base_dir + "/figures")

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Time to 1 Significant Quantity

fig, ax = plt.subplots()
ax.spines["top"].set_color("None")
ax.spines["right"].set_color("None")

ax.scatter(masses, U_time_to_SQ/24, label="$^{238}$U", marker='o', color='r')
ax.scatter(masses, Th_time_to_SQ/24, label="$^{232}$Th", marker='s', color='g')

ax.set_ylim(10, 200)

np.save("U_time_to_SQ_depletion", U_time_to_SQ)
np.save("Th_time_to_SQ_depletion", Th_time_to_SQ)

# Fit data to 1/x function:
def fit(x, A, B):
    return (A/x) + B

U_popt, U_pcov = curve_fit(fit, masses, U_time_to_SQ)
Th_popt, Th_pcov = curve_fit(fit, masses, Th_time_to_SQ)

fit_masses = np.linspace(1, masses[-1], num=100)
ax.plot(fit_masses, fit(fit_masses, *U_popt)/24, alpha=0.3, color='r')
ax.plot(fit_masses, fit(fit_masses, *Th_popt)/24, alpha=0.3, color='g')

ax.legend()

ax.set_xlim(0, masses[-1] + 2)
ax.set_ylim(10, np.max(Th_time_to_SQ/24) + 100)

ax.set_yscale("log")

ax.set_title("Time to Breed a Significant Quantity of Fissile Material", fontsize=14)
ax.set_ylabel("Time (days)", fontsize=14)
ax.set_xlabel("Mass of Fertile Material (metric tons)", fontsize=14)

fig.savefig("time_to_sq.png", dpi=300)

for i, mass in enumerate(masses):
    fig, ax = plt.subplots()

    ax.plot(time_steps, U_fissile_masses[i], label="Pu239")
    ax.plot(time_steps, Th_fissile_masses[i], label="U233")

    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Mass (kg)")
    ax.set_title("Fissile Mass vs. Time for a Fertile Mass of " + str(mass) + " metric tons")

    ax.legend()

    fig.savefig(str(mass) + "_metric_tons.png", dpi=300)
