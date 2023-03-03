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
print(masses)
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

# ====================================================
# Fission Power
# ====================================================
U_fission_powers = np.empty((len(masses), num_steps, 2))
Th_fission_powers = np.empty((len(masses), num_steps, 2))

for i, mass in enumerate(masses):
    """ Uranium """
    os.chdir(base_dir + "/Uranium/" + str(mass))

    for step in range(0, num_steps):
        sp = openmc.StatePoint('openmc_simulation_n'+str(step)+'.h5')
        tally = sp.get_tally(name='FLiBe Tally')

        U_fission_powers[i, step, 0] = anp.get_uvalue(tally, 'kappa-fission').n * total_neutron_rate * MJ_per_eV
        U_fission_powers[i, step, 1] = anp.get_uvalue(tally, 'kappa-fission').s * total_neutron_rate * MJ_per_eV

    os.chdir('../../..')

    """ Thorium """
    os.chdir(base_dir + "/Thorium/" + str(mass))

    for step in range(0, num_steps):
        sp = openmc.StatePoint('openmc_simulation_n'+str(step)+'.h5')
        tally = sp.get_tally(name='FLiBe Tally')

        Th_fission_powers[i, step, 0] = anp.get_uvalue(tally, 'kappa-fission').n * total_neutron_rate * MJ_per_eV
        Th_fission_powers[i, step, 1] = anp.get_uvalue(tally, 'kappa-fission').s * total_neutron_rate * MJ_per_eV

    os.chdir('../../..')



# ====================================================
# Decay Heat
# ====================================================


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

""" Time to 1 Significant Quantity """
fig, ax = plt.subplots()
ax.spines["top"].set_color("None")
ax.spines["right"].set_color("None")

ax.scatter(masses, U_time_to_SQ/24, label="$^{238}$U", marker='o')
ax.scatter(masses, Th_time_to_SQ/24, label="$^{232}$Th", marker='s')

np.save("U_time_to_SQ_depletion", U_time_to_SQ)
np.save("Th_time_to_SQ_depletion", Th_time_to_SQ)

# Fit data to 1/x function:
def fit(x, A, B):
    return A/x + B

U_popt, U_pcov = curve_fit(fit, masses, U_time_to_SQ)
Th_popt, Th_pcov = curve_fit(fit, masses, Th_time_to_SQ)

fit_masses = np.linspace(1, masses[-1], num=100)
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

# Fission Power:
fig, ax = plt.subplots()
ax.spines["top"].set_color("None")
ax.spines["right"].set_color("None")

U_fission_power_at_SQ = np.empty((len(masses)))
Th_fission_power_at_SQ = np.empty((len(masses)))
for i, mass in enumerate(masses):
    """ Uranium """
    U_res = linregress(time_steps, U_fission_powers[i, :, 0])
    U_fission_power_at_SQ[i] = U_res.intercept + U_res.slope*U_time_to_SQ[i]

    """ Thorium """
    Th_res = linregress(time_steps, Th_fission_powers[i, :, 0])
    Th_fission_power_at_SQ[i] = Th_res.intercept + Th_res.slope*Th_time_to_SQ[i]

ax.scatter(masses, U_fission_powers[:, 0, 0], c='r', s=4, label='At t = 0')
ax.scatter(masses, Th_fission_powers[:, 0, 0], c='g', s=4, label='At t = 0')

ax.scatter(masses, U_fission_power_at_SQ, c='r', s=4, label='After 1 SQ bred')
ax.scatter(masses, Th_fission_power_at_SQ, c='g', s=4, label='After 1 SQ bred')

ax.fill_between(masses, U_fission_powers[:, 0, 0], U_fission_power_at_SQ, color='r', alpha=0.3)
ax.fill_between(masses, Th_fission_powers[:, 0, 0], Th_fission_power_at_SQ, color='g', alpha=0.3)

text_offset = 5
ax.annotate("t = $t_{SQ}$", (masses[-1], U_fission_power_at_SQ[-1]), color='r', textcoords='offset points', xytext=(text_offset, 0))
ax.annotate("t = 0", (masses[-1], U_fission_powers[-1, 0, 0]), color='r', textcoords='offset points', xytext=(text_offset, 0))

ax.annotate("t = $t_{SQ}$", (masses[-1], Th_fission_power_at_SQ[-1]), color='g', textcoords='offset points', xytext=(text_offset, 0))
ax.annotate("t = 0", (masses[-1], Th_fission_powers[-1, 0, 0]), color='g', textcoords='offset points', xytext=(text_offset, 0))

ax.set_xlim(0, masses[-1] + 5)
ax.set_ylim(0, U_fission_power_at_SQ[-1] + 10)

ax.set_title("Fission Power in Doped FLiBe Blanket", fontsize=14)
ax.set_ylabel("Fission Power (MW)", fontsize=14)
ax.set_xlabel("Fertile Mass (metric tons)", fontsize=14)

fig.savefig("fission_power.png")