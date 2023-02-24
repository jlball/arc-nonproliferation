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

    breeding_rate = U_tally.get_values(scores=['absorption'], nuclides=['U238']) * total_neutron_rate * Pu239_mass_in_kg
    U_time_to_SQ[i] = sig_quantity / breeding_rate # seconds
    os.chdir("../../..")

    """ Extract time to 1 SQ for Thorium """
    os.chdir(base_dir + "/Thorium/" + str(mass))
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

U_tbr = np.empty(len(masses))
Th_tbr = np.empty(len(masses))

for i, mass in enumerate(masses):

    """ Extract time to 1 SQ for Uranium """
    os.chdir(base_dir + "/Uranium/" + str(mass))
    sp = openmc.StatePoint('statepoint.10.h5')
    U_tally = sp.get_tally(name='FLiBe Tally')

    U_tbr[i] = U_tally.get_values(scores=['(n,Xt)'])
    os.chdir("../../..")

    """ Extract time to 1 SQ for Thorium """
    os.chdir(base_dir + "/Thorium/" + str(mass))
    sp = openmc.StatePoint('statepoint.10.h5')
    Th_tally = sp.get_tally(name='FLiBe Tally')

    Th_tbr[i] = Th_tally.get_values(scores=['(n,Xt)'])
    os.chdir("../../..")

# ====================================================
# Plotting
# ====================================================

masses = masses/1e3

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

ax.scatter(masses, U_time_to_SQ/24, label="$^{238}$U", marker='o')
ax.scatter(masses, Th_time_to_SQ/24, label="$^{232}$Th", marker='s')

# Fit data to 1/x function:
def fit(x, A, B):
    return A/x + B

U_popt, U_pcov = curve_fit(fit, masses, U_time_to_SQ)
Th_popt, Th_pcov = curve_fit(fit, masses, Th_time_to_SQ)

fit_masses = np.linspace(5, masses[-1], num=100)
ax.plot(fit_masses, fit(fit_masses, *U_popt)/24, alpha=0.3)
ax.plot(fit_masses, fit(fit_masses, *Th_popt)/24, alpha=0.3)

ax.legend()

ax.set_xlim(0, masses[-1] + 2)
ax.set_ylim(0, np.max(Th_time_to_SQ/24) + 5)

ax.set_title("Time to Breed a Significant Quantity of Fissile Material", fontsize=14)
ax.set_ylabel("Time (days)", fontsize=14)
ax.set_xlabel("Mass of Fertile Material (metric tons)", fontsize=14)

fig.savefig("time_to_sq.png")

# Fissile Proliferance

fusion_power = 500 #MW

U_fissile_proliferance = 1/(U_time_to_SQ * masses * fusion_power)
Th_fissile_proliferance = 1/(Th_time_to_SQ * masses * fusion_power)

fig, ax = plt.subplots()

ax.scatter(masses, U_fissile_proliferance)
ax.scatter(masses, Th_fissile_proliferance)

fig.savefig("fissile_proliferance.png")

# TBR
fig, ax = plt.subplots()
ax.spines["top"].set_color("None")
ax.spines["right"].set_color("None")

ax.scatter(masses, U_tbr, label="$^{238}$U", marker='o')
ax.scatter(masses, Th_tbr, label="$^{232}$Th", marker='s')

U_res = linregress(masses, y=U_tbr)
Th_res = linregress(masses, y=Th_tbr)

ax.plot(masses, U_res.intercept + U_res.slope*masses, alpha=0.3)
ax.plot(masses, Th_res.intercept + Th_res.slope*masses, alpha=0.3)

ax.legend()

ax.set_title("Tritium Breeding Ratio", fontsize=14)
ax.set_ylabel("TBR", fontsize=14)
ax.set_xlabel("Mass of Fertile Material (metric tons)", fontsize=14)

ax.set_xlim(0, masses[-1] + 2)
ax.set_ylim(0.6, 1.2)

fig.savefig("tbr.png")