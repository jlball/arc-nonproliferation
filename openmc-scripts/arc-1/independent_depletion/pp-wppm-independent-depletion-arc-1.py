import openmc
from openmc.deplete import Results
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mpl
import numpy as np
import sys
from arc_nonproliferation.postprocess import *
from arc_nonproliferation.constants import *
from scipy.optimize import curve_fit
from scipy.stats import linregress
from numpy.polynomial.polynomial import Polynomial
import time

"""
This script generates figures from data generated by the script:

wppm-independent-depletion-arc-1.py

This script specifically focuses on the impact of uranium impurities
in FLiBe.

Usage: Call the script from the command line, with the additional
command line argument of the name of the output directory produced
by independent-depletion-arc-1.py you wish to analyze:

python3 pp-independent-depletion-arc-1.py <output_directory_name>

"""

# get base directory from command line argument
if sys.argv[1] is not None:
    base_dir = './' + sys.argv[1]
else:
    raise ValueError("No base directory specified!")

openmc.config['chain_file'] = anp.constants.chain_file

# ====================================================
# Extract Data
# ====================================================

""" Load masses and initialisze final output arrays """
masses = np.loadtxt(base_dir + '/masses.txt')

os.chdir(base_dir + "/Uranium/" + str(int(masses[0])))
time_steps = Results("depletion_results.h5").get_times()
os.chdir("../../..")

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Fissile Mass

""" Iterate through each mass simulated and get fissile mass at each time step"""
U_fissile_masses = np.empty((len(masses), len(time_steps)))
U_fertile_masses = np.empty((len(masses), len(time_steps)))

for i, mass in enumerate(masses):

    os.chdir(base_dir + "/Uranium/" + str(int(mass)))

    U_results = Results('depletion_results.h5')
    
    #While we're here, get the number of depletion steps:
    time_steps = U_results.get_times()
    num_steps = len(time_steps)

    U_fissile_masses[i] = get_masses_from_mats('Pu239', U_results, density=True)*1e6 # convert from kg/cm3 to kg/m3
    U_fertile_masses[i] = get_masses_from_mats('U238', U_results, density=True)*1e6 # convert from kg/cm3 to kg/m3

    mats = openmc.Materials.from_xml()
    #wppm = (mats[2].get_mass("U238")/mats[2].get_mass())*1e6
    #print(f"specified wppm: {mass} true wppm: {wppm}")
    U_mass = 0
    U_mass += get_element_mass(mats[2], "U")

    print(f"Uranium mass:{U_mass/1e3} (kg) for impurity level:{mass}")

    os.chdir("../../..")

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Decay Photon Spectrum

# U_decay_spectra_channels = []
# U_decay_spectra_blanket = []

# Th_decay_spectra_channels = []
# Th_decay_spectra_blanket = []

# for i in range(0, num_steps):
#     """ Uranium """
#     os.chdir(base_dir + "/Uranium/" + str(mass))

#     U_results = Results('depletion_results.h5')
#     U_mats = U_results.export_to_materials(i)

#     flibe_mat_channels = get_material_by_name(U_mats, "doped flibe channels")
#     flibe_mat_blanket = get_material_by_name(U_mats, "doped flibe blanket")

#     U_decay_spectra_channels.append(flibe_mat_channels.decay_photon_energy)
#     U_decay_spectra_blanket.append(flibe_mat_blanket.decay_photon_energy)

#     os.chdir('../../..')

#     """ Thorium """
#     os.chdir(base_dir + "/Thorium/" + str(mass))

#     Th_results = Results('depletion_results.h5')
#     Th_mats = Th_results.export_to_materials(i)

#     flibe_mat_channels = get_material_by_name(Th_mats, "doped flibe channels")
#     flibe_mat_blanket = get_material_by_name(Th_mats, "doped flibe blanket")

#     Th_decay_spectra_channels.append(flibe_mat_channels.decay_photon_energy)
#     Th_decay_spectra_blanket.append(flibe_mat_blanket.decay_photon_energy)

#     os.chdir('../../..')

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Isotopic Purity

# U_purities = np.empty(len(masses))
# Th_purities = np.empty(len(masses))

# for i, mass in enumerate(masses):
#     """ Uranium """
#     os.chdir(base_dir + "/Uranium/" + str(mass))

#     U_results = Results('depletion_results.h5')
#     U_purity = extract_isotopic_purity("U", U_results)
#     U_purities[i] = U_purity[-1]

#     os.chdir('../../..')

#     """ Thorium """
#     os.chdir(base_dir + "/Thorium/" + str(mass))

#     Th_results = Results('depletion_results.h5')
#     Th_purity = extract_isotopic_purity("Th", Th_results)
#     Th_purities[i] = Th_purity[-1]

#     os.chdir('../../..')

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Contact dose rate
    
# dose_rate_steps = 25
# contact_dose_rate = np.empty((len(masses), dose_rate_steps))
# dose_rate_times = np.empty(dose_rate_steps)
# for i, mass in enumerate(masses):

#     os.chdir(base_dir + "/Uranium/" + str(int(mass)))

#     U_results = Results('depletion_results.h5')

#     for j in range(0, dose_rate_steps):
#         start_time = time.perf_counter()

#         mats = U_results.export_to_materials(j*2)
#         dose_rate_times[j] = time_steps[j*2]

#         channel_mat = get_material_by_name(mats, "doped flibe channels")
#         blanket_mat = get_material_by_name(mats, "doped flibe blanket")

#         contact_dose_rate[i, j] = extract_contact_dose_rate(blanket_mat) + extract_contact_dose_rate(channel_mat)
#         print(f"dose rate function time: {time.perf_counter() - start_time}")

#     os.chdir("../../..")
# ====================================================
# Plotting
# ====================================================

width_in = 7
height_in = 5

fontdict = {"size":16}

tick_font_size = 12

#Change into dedicated directory for figures or create figures directory
try:
    os.chdir(base_dir + "/figures")
except:
    os.mkdir(base_dir + "/figures")
    os.chdir(base_dir + "/figures")

dpi = 300

title_y = 1.05

fontdict = {"size":16}

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Fissile Mass

fig_tot, ax_tot = plt.subplots()
ax_tot.spines["top"].set_color("None")
ax_tot.spines["right"].set_color("None")

fig_tot.set_size_inches(width_in, height_in)

U_cm = mpl.cm.Reds
Th_cm = mpl.cm.Greens

norm = colors.Normalize(vmin=-10, vmax=masses.max() + 10)

for i, mass in enumerate(masses):
    # Individual plot for this mass
    fig, ax = plt.subplots()

    ax.scatter(time_steps/365, U_fissile_masses[i]*1e3, label="Pu239")
    ax.scatter(time_steps/365, U_fertile_masses[i]*1e3, label="U238")

    ax.set_xlabel("Time (years)", fontdict=fontdict)
    ax.set_ylabel("g/m$^3$", fontdict=fontdict)
    ax.set_title(f"Fissile Mass vs. Time for a beryllium impurity concentration of {mass} wppm")

    ax.legend()

    fig.savefig(str(mass) + "_metric_tons.png", dpi=300)

    # Add this mass to the total plot
    U_color = U_cm.__call__(norm(mass))
    Th_color = Th_cm.__call__(norm(mass))

    ax_tot.plot(time_steps/365, U_fissile_masses[i]*1e3, label=str(mass), color=U_color)
    #ax_tot.plot(time_steps/365, U_fertile_masses[i]*1e3, color=U_color)
    ax_tot.annotate(f"{mass} wppm", (time_steps[-1]/365, U_fissile_masses[i, -1]*1e3), color=U_color, textcoords="offset points", xytext=(-60, 10))

ax_tot.set_xlabel("Time (years)", fontdict=fontdict)
ax_tot.set_ylabel("g/m$^3$", fontdict=fontdict)
ax_tot.set_title("Mass of Pu-239 vs. Time", fontdict=fontdict, y=title_y)

ax_tot.set_ylim(0, np.max(1.05*U_fissile_masses*1e3))
ax_tot.set_xlim(0, time_steps[-1]/365)

ax_tot.set_yscale("linear")
#ax_tot.tight_layout()

fig_tot.savefig("all_masses.png", dpi=dpi)

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Fertile / Fissile Ratio

fig, ax = plt.subplots()
ax.spines["top"].set_color("None")
ax.spines["right"].set_color("None")

fig.set_size_inches(width_in, height_in)

for i, mass in enumerate(masses):
    U_color = U_cm.__call__(norm(mass))
    ax.plot(time_steps/365, U_fertile_masses[i]/U_fissile_masses[i], color=U_color)
    ax.annotate(f"{mass} wppm", (time_steps[-1]/365, U_fertile_masses[i, -1]/U_fissile_masses[i, -1]), color=U_color, textcoords="offset points", xytext=(5, 0))

ax.set_xlabel("Time (years)", fontdict=fontdict)
#ax_tot.set_ylabel("", fontdict=fontdict)
ax.set_title("Mass ratio of U-238 to Pu-239", fontdict=fontdict, y=title_y)

fig.tight_layout()
fig.savefig("fertile_fissile_ratio.png", dpi=dpi)

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Contact dose rate

# fig, ax = plt.subplots()
# ax.spines["top"].set_color("None")
# ax.spines["right"].set_color("None")

# fig.set_size_inches(width_in, height_in)

# for i, mass in enumerate(masses):
#     U_color = U_cm.__call__(norm(mass))
#     ax.plot(dose_rate_times/365, contact_dose_rate[i], label=mass, color=U_color)
#     ax.annotate(f"{mass} wppm", (dose_rate_times[-1]/365, contact_dose_rate[i, -1]), color=U_color, textcoords="offset points", xytext=(5, 0))

#     ax.set_xlabel("Time (years)", fontdict=fontdict)
#     ax.set_ylabel("Contact dose rate (Sv/hr)", fontdict=fontdict)
#     ax.set_title("Contact dose rate in impure FLiBe", fontdict=fontdict)

# fig.tight_layout()
# fig.savefig("contact_dose_rate.png", dpi=300)

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Decay Photon Spectrum

# Uranium

# fig, axs = plt.subplots(1, 2)

# for ax in axs:
#     ax.spines["top"].set_color("None")
#     ax.spines["right"].set_color("None")

#     ax.set_yscale("log")

#     ax.set_ylim(1e15, 1e19)
#     ax.set_xlim(0, 3)
#     ax.set_xlabel("Photon Energy (MeV)")

# for i in range(0, num_steps):
#     axs[0].step(U_decay_spectra_channels[i].x/1e6, U_decay_spectra_channels[i].p, label = str("Step " + str(int(i))))
#     axs[1].step(U_decay_spectra_blanket[i].x/1e6, U_decay_spectra_blanket[i].p, label = str("Step " + str(int(i))))
    
# axs[0].set_title("Gamma spectrum in a Uranium doped cooling channel")
# axs[1].set_title("Gamma spectrum in a Uranium doped blanket")

# fig.set_size_inches(10, 4)
# fig.savefig("U_decay_spectra.png")

# Thorium:

# fig, ax = plt.subplots()
# ax.spines["top"].set_color("None")
# ax.spines["right"].set_color("None")

# for i, dist in enumerate(Th_decay_spectra):
#     ax.step(dist.x, dist.p, label = str("Step " + str(int(i))))

# ax.set_yscale("log")

# ax.set_title("Gamma spectrum at $t_{SQ}$ in a Thorium doped blanket")
# ax.set_xlabel("Photon Energy (eV)")

# ax.set_ylim(1e15, 1e22)

# fig.savefig("Th_decay_spectra.png")

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Isotopic Purity

# fig, ax = plt.subplots()
# ax.spines["top"].set_color("None")
# ax.spines["right"].set_color("None")

# ax.scatter(masses, U_purities*100, label = "Pu239", color='r')
# ax.scatter(masses, Th_purities*100, label = "U233", color='g')

# ax.legend()

# ax.set_title("Isotopic Purity vs. Fertile Inventory", fontsize=14)
# ax.set_ylabel("Isotopic Purity (\% fissile isotope)", fontsize=14)
# ax.set_xlabel("Fertile Mass (metric tons)", fontsize=14)

# fig.savefig("isotopic_purity.png")


