import openmc
from openmc.deplete import Results
import matplotlib.pyplot as plt
import numpy as np
import sys
from arc_nonproliferation.postprocess import *
from scipy.optimize import curve_fit

if sys.argv[1] is not None:
    base_dir = './' + sys.argv[1]
else:
    raise ValueError("No base directory specified!")

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
    print("============== URANIUM MASS:" + str(mass) + "kg ==============")
    os.chdir(base_dir + "/Uranium/" + str(mass))

    U_results = Results('depletion_results.h5')
    U_time_to_SQ[i] = extract_time_to_sq('U', U_results)

    os.chdir("../../..")

    """ Extract time to 1 SQ for Thorium """
    print("============== THORIUM MASS:" + str(mass) + "kg ==============")
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

# Fit data to 1/x function:
def fit(x, A, B):
    return A/x + B

U_popt, U_pcov = curve_fit(fit, masses, U_time_to_SQ)
Th_popt, Th_pcov = curve_fit(fit, masses, Th_time_to_SQ)

fit_masses = np.linspace(5, masses[-1], num=100)
ax.plot(fit_masses/1e3, fit(fit_masses, *U_popt)/24, alpha=0.3)
ax.plot(fit_masses/1e3, fit(fit_masses, *Th_popt)/24, alpha=0.3)

ax.legend()

ax.set_xlim(0, masses[-1]/1e3 + 2)
ax.set_ylim(0, np.max(Th_time_to_SQ/24) + 5)

ax.set_title("Time to Breed a Significant Quantity of Fissile Material", fontsize=14)
ax.set_ylabel("Time (days)", fontsize=14)
ax.set_xlabel("Mass of Fertile Material (metric tons)", fontsize=14)

fig.savefig("time_to_sq.png")

# Fissile Proliferance

fusion_power = 500 #MW

U_fissile_proliferance = 1/(U_time_to_SQ * masses/1e3 * fusion_power)
Th_fissile_proliferance = 1/(Th_time_to_SQ * masses/1e3 * fusion_power)

fig, ax = plt.subplots()

ax.scatter(masses, U_fissile_proliferance)
ax.scatter(masses, Th_fissile_proliferance)

fig.savefig("fissile_proliferance.png")