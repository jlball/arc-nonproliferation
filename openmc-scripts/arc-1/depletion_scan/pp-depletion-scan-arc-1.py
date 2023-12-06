import openmc
from openmc.deplete import Results
import matplotlib.pyplot as plt
import numpy as np
import sys
from arc_nonproliferation.postprocess import *
from arc_nonproliferation.constants import *
from scipy.optimize import curve_fit
from scipy.stats import linregress
from scipy.interpolate import interp1d
import time
import pickle

if sys.argv[1] is not None:
    base_dir = './' + sys.argv[1]
else:
    raise ValueError("No base directory specified!")

fusion_power = 500 #MW
total_neutron_rate = fusion_power * neutrons_per_MJ

openmc.config['chain_file'] = chain_file

# ====================================================
# Extract Data
# ====================================================

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Time to a Significant Quantity

init_time = time.perf_counter()

""" Load masses and initialisze final output arrays """
masses = np.loadtxt(base_dir + '/masses.txt')

U_time_to_SQ = np.empty(len(masses))
Th_time_to_SQ = np.empty(len(masses))

U_tsq_idx = np.empty(len(masses))
Th_tsq_idx = np.empty(len(masses))

""" Iterate through each mass simulated and compute time to SQ"""
for i, mass in enumerate(masses):

    """ Extract time to 1 SQ for Uranium """
    os.chdir(base_dir + "/Uranium/" + str(mass))

    U_results = Results('depletion_results.h5')
    U_tsq_idx[i], U_time_to_SQ[i] = extract_time_to_sq('U', U_results)

    #While we're here, get the number of depletion steps:
    U_time_steps = U_results.get_times()
    num_steps = len(U_time_steps)

    os.chdir("../../..")

    """ Extract time to 1 SQ for Thorium """
    os.chdir(base_dir + "/Thorium/" + str(mass))

    Th_results = Results('depletion_results.h5')
    Th_tsq_idx[i], Th_time_to_SQ[i] = extract_time_to_sq('Th', Th_results)

    Th_time_steps = Th_results.get_times()

    os.chdir("../../..")

print("Loaded time to 1 SQ data in " + str(round(time.perf_counter() - init_time, 2)) + " seconds.")

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Fission
init_time = time.perf_counter()

U_fission_powers = np.empty((len(masses), num_steps))
Th_fission_powers = np.empty((len(masses), num_steps))

U_fission_rates = np.empty((len(masses), num_steps))
Th_fission_rates = np.empty((len(masses), num_steps))

for i, mass in enumerate(masses):
    """ Uranium """
    os.chdir(base_dir + "/Uranium/" + str(mass))

    for step in range(0, num_steps):
        sp = openmc.StatePoint('openmc_simulation_n'+str(step)+'.h5')
        tally = sp.get_tally(name='FLiBe Tally')

        U_fission_powers[i, step] = np.sum(tally.get_values(scores=["kappa-fission"])) * total_neutron_rate * MJ_per_eV
        
    os.chdir('../../..')

    """ Thorium """
    os.chdir(base_dir + "/Thorium/" + str(mass))

    for step in range(0, num_steps):
        sp = openmc.StatePoint('openmc_simulation_n'+str(step)+'.h5')
        tally = sp.get_tally(name='FLiBe Tally')

        Th_fission_powers[i, step] = np.sum(tally.get_values(scores=["kappa-fission"])) * total_neutron_rate * MJ_per_eV

    os.chdir('../../..')

print("Loaded fission power data in " + str(round(time.perf_counter() - init_time, 2)) + " seconds.")

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# TBR
init_time = time.perf_counter()

U_TBR = np.empty((len(masses), 2, 2))
Th_TBR = np.empty((len(masses), 2, 2))

for i, mass in enumerate(masses):
    """ Uranium """
    os.chdir(base_dir + "/Uranium/" + str(mass))

    sp = openmc.StatePoint('openmc_simulation_n0.h5')
    tally = sp.get_tally(name='FLiBe Tally')

    U_TBR[i, 0, 0] = np.sum(tally.get_values(scores=["(n,Xt)"]))
    U_TBR[i, 0, 1] = np.sqrt(np.sum(np.square(tally.get_values(scores=["(n,Xt)"], value="std_dev"))))
    sp.close()

    sp = openmc.StatePoint('openmc_simulation_n' + str(num_steps-1) + '.h5')
    tally = sp.get_tally(name='FLiBe Tally')

    U_TBR[i, 1, 0] = np.sum(tally.get_values(scores=["(n,Xt)"]))
    U_TBR[i, 1, 1] = np.sqrt(np.sum(np.square(tally.get_values(scores=["(n,Xt)"], value="std_dev"))))
    sp.close()
        
    os.chdir('../../..')

    """ Thorium """
    os.chdir(base_dir + "/Thorium/" + str(mass))

    sp = openmc.StatePoint('openmc_simulation_n0.h5')
    tally = sp.get_tally(name='FLiBe Tally')

    Th_TBR[i, 0, 0] = np.sum(tally.get_values(scores=["(n,Xt)"]))
    Th_TBR[i, 0, 1] = np.sqrt(np.sum(np.square(tally.get_values(scores=["(n,Xt)"], value="std_dev"))))
    sp.close()

    sp = openmc.StatePoint('openmc_simulation_n' + str(num_steps-1) + '.h5')
    tally = sp.get_tally(name='FLiBe Tally')

    Th_TBR[i, 1, 0] = np.sum(tally.get_values(scores=["(n,Xt)"]))
    Th_TBR[i, 1, 1] = np.sqrt(np.sum(np.square(tally.get_values(scores=["(n,Xt)"], value="std_dev"))))
    sp.close()

    os.chdir('../../..')

print("Loaded TBR data in " + str(round(time.perf_counter() - init_time, 2)) + " seconds.")

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Flux Spectrum
init_time = time.perf_counter()

U_flux_spectra = np.empty((len(masses), num_steps, 709, 2))
Th_flux_spectra = np.empty((len(masses), num_steps, 709, 2))

""" Uranium """
for i, mass in enumerate(masses):
    """ Uranium """
    os.chdir(base_dir + "/Uranium/" + str(mass))

    for step in range(0, num_steps):
        sp = openmc.StatePoint('openmc_simulation_n'+str(step)+'.h5')
        flux_tally = sp.get_tally(name='Channel Flux Tally')
        flux_spectrum = flux_tally.get_reshaped_data()
        U_flux_spectra[i, step] = flux_spectrum.reshape((709,2))

    os.chdir('../../..')

    """ Thorium """
    os.chdir(base_dir + "/Thorium/" + str(mass))

    for step in range(0, num_steps):
        sp = openmc.StatePoint('openmc_simulation_n'+str(step)+'.h5')
        flux_tally = sp.get_tally(name='Channel Flux Tally')
        flux_spectrum = flux_tally.get_reshaped_data()
        Th_flux_spectra[i, step] = flux_spectrum.reshape((709,2))

    os.chdir('../../..')

energy_groups = openmc.mgxs.EnergyGroups(openmc.mgxs.GROUP_STRUCTURES['CCFE-709'])
energies = energy_groups.group_edges

print("Loaded flux spectrum data in "  + str(round(time.perf_counter() - init_time, 2)) + " seconds.")

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Isotopic Purity
init_time = time.perf_counter()

U_purities = np.empty(len(masses))
Th_purities = np.empty(len(masses))

#fig, ax = plt.subplots()
for i, mass in enumerate(masses):
    """ Uranium """
    os.chdir(base_dir + "/Uranium/" + str(mass))

    U_results = Results('depletion_results.h5')
    U_purity = extract_isotopic_purity("U", U_results)

    #ax.plot(U_time_steps[1:], U_purity[1:])
    #ax.vlines(U_time_to_SQ[i]/24, U_purity[1:].min(), U_purity[1:].max())

    # Use linear fit to extract purity at t_SQ
    U_purity_fit = linregress(U_time_steps[1:], y=U_purity[1:])
    U_purities[i] = U_purity_fit.slope * (U_time_to_SQ[i]/24) + U_purity_fit.intercept

    os.chdir('../../..')

    """ Thorium """
    os.chdir(base_dir + "/Thorium/" + str(mass))

    Th_results = Results('depletion_results.h5')
    Th_purity = extract_isotopic_purity("Th", Th_results)

    # Use linear fit to extract purity at t_SQ
    Th_purity_fit = linregress(Th_time_steps[1:], y=Th_purity[1:])
    Th_purities[i] = Th_purity_fit.slope * (Th_time_to_SQ[i]/24) + Th_purity_fit.intercept

    os.chdir('../../..')

print("Loaded isotopic purity data in "  + str(round(time.perf_counter() - init_time, 2)) + "seconds.")

#ax.scatter(U_time_to_SQ/24, U_purities)

#plt.show()

# # +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Decay Photon Spectrum
init_time = time.perf_counter()

U_decay_spectra_channels = []
U_decay_spectra_blanket = []

Th_decay_spectra_channels = []
Th_decay_spectra_blanket = []

for i, mass in enumerate(masses):
    """ Uranium """
    os.chdir(base_dir + "/Uranium/" + str(mass))

    U_results = Results('depletion_results.h5')
    U_mats = U_results.export_to_materials(int(U_tsq_idx[i]))

    # Get materials at time step just before and just after t_SQ
    mat_channels = get_material_by_name(U_mats, "doped flibe channels")
    mat_blanket = get_material_by_name(U_mats, "doped flibe blanket")

    U_decay_spectra_channels.append(mat_channels.get_decay_photon_energy(units="Bq/cm3"))
    U_decay_spectra_blanket.append(mat_blanket.get_decay_photon_energy(units="Bq/cm3"))

    os.chdir('../../..')

    """ Thorium """
    os.chdir(base_dir + "/Thorium/" + str(mass))

    Th_results = Results('depletion_results.h5')
    Th_mats = Th_results.export_to_materials(int(Th_tsq_idx[i]))

    flibe_mat_channels = get_material_by_name(Th_mats, "doped flibe channels")
    flibe_mat_blanket = get_material_by_name(Th_mats, "doped flibe blanket")

    Th_decay_spectra_channels.append(flibe_mat_channels.get_decay_photon_energy(units="Bq/cm3"))
    Th_decay_spectra_blanket.append(flibe_mat_blanket.get_decay_photon_energy(units="Bq/cm3"))

    os.chdir('../../..')

print("Loaded decay photon data in "  + str(round(time.perf_counter() - init_time, 2)) + "seconds.")

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Fissile Mass
init_time = time.perf_counter()

""" Iterate through each mass simulated and get fissile mass at each time step"""
U_fissile_masses = np.empty((len(masses), len(U_time_steps)))
Th_fissile_masses = np.empty((len(masses), len(Th_time_steps)))
for i, mass in enumerate(masses):

    os.chdir(base_dir + "/Uranium/" + str(mass))

    U_results = Results('depletion_results.h5')
    U_fissile_masses[i] = get_masses_from_mats('Pu239', U_results)

    os.chdir("../../..")

    os.chdir(base_dir + "/Thorium/" + str(mass))

    Th_results = Results('depletion_results.h5')
    Th_fissile_masses[i] = get_masses_from_mats('U233', Th_results)

    os.chdir("../../..")

print("Loaded fissile mass data in "  + str(round(time.perf_counter() - init_time, 2)) + " seconds.")

# # +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Reaction Spectra (energy resolved n,gamma and fission on fertile isotope)

init_time = time.perf_counter()

U_reaction_spectra = np.empty((len(masses), num_steps, 709, 2, 2))
Th_reaction_spectra = np.empty((len(masses), num_steps, 709, 2, 2))

""" Uranium """
for i, mass in enumerate(masses):
    os.chdir(base_dir + "/Uranium/" + str(mass))

    for step in range(0, num_steps):
        sp = openmc.StatePoint('openmc_simulation_n'+str(step)+'.h5')
        reaction_spectra_tally = sp.get_tally(name='Fertile Tally')
        reaction_spectra = reaction_spectra_tally.get_reshaped_data()
        U_reaction_spectra[i, step] = reaction_spectra.reshape((709,2,2))

    os.chdir('../../..')

""" Thorium """
for i, mass in enumerate(masses):
    os.chdir(base_dir + "/Thorium/" + str(mass))

    for step in range(0, num_steps):
        sp = openmc.StatePoint('openmc_simulation_n'+str(step)+'.h5')
        absorption_tally = sp.get_tally(name='Fertile Tally')
        absorption_spectrum = absorption_tally.get_reshaped_data()
        Th_reaction_spectra[i, step] = absorption_spectrum.reshape((709,2,2))

    os.chdir('../../..')

print("Loaded reaction spectra data in "  + str(round(time.perf_counter() - init_time, 2)) + " seconds.")

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Decay Heat
init_time = time.perf_counter()

U_decay_heats = np.empty(len(masses))
Th_decay_heats = np.empty(len(masses))

""" Uranium """
for i, mass in enumerate(masses):
    os.chdir(base_dir + "/Uranium/" + str(mass))

    U_results = Results('depletion_results.h5')

    idx = np.abs(U_fissile_masses[i] - anp.sig_quantity).argmin()

    if U_fissile_masses[i, idx] < anp.sig_quantity:
        U_mats_pre_tsq = U_results.export_to_materials(idx)
        U_mats_post_tsq = U_results.export_to_materials(idx+1)
    else:
        U_mats_pre_tsq = U_results.export_to_materials(idx-1)
        U_mats_post_tsq = U_results.export_to_materials(idx)

    # Get materials at time step just before and just after t_SQ
    mat_channels_pre_tsq = get_material_by_name(U_mats_pre_tsq, "doped flibe channels")
    mat_blanket_pre_tsq = get_material_by_name(U_mats_pre_tsq, "doped flibe blanket")

    mat_channels_post_tsq = get_material_by_name(U_mats_post_tsq, "doped flibe channels")
    mat_blanket_post_tsq = get_material_by_name(U_mats_post_tsq, "doped flibe blanket")

    decay_heat_pre_tsq =  mat_channels_pre_tsq.get_decay_heat() + mat_blanket_pre_tsq.get_decay_heat()
    decay_heat_post_tsq =  mat_channels_post_tsq.get_decay_heat() + mat_blanket_post_tsq.get_decay_heat()

    if U_fissile_masses[i, idx] < anp.sig_quantity:
        U_decay_heats[i] = np.interp(U_time_to_SQ[i], U_time_steps[idx:idx+2], [decay_heat_pre_tsq, decay_heat_post_tsq])
    else:
        U_decay_heats[i] = np.interp(U_time_to_SQ[i], U_time_steps[idx-1:idx+1], [decay_heat_pre_tsq, decay_heat_post_tsq])

    os.chdir('../../..')

""" Thorium """
for i, mass in enumerate(masses):
    os.chdir(base_dir + "/Thorium/" + str(mass))

    Th_results = Results('depletion_results.h5')
    idx = np.abs(Th_fissile_masses[i] - anp.sig_quantity).argmin()

    if Th_fissile_masses[i, idx] < anp.sig_quantity:
        Th_mats_pre_tsq = Th_results.export_to_materials(idx)
        Th_mats_post_tsq = Th_results.export_to_materials(idx+1)
    else:
        Th_mats_pre_tsq = Th_results.export_to_materials(idx-1)
        Th_mats_post_tsq = Th_results.export_to_materials(idx)

    # Get materials at time step just before and just after t_SQ
    mat_channels_pre_tsq = get_material_by_name(Th_mats_pre_tsq, "doped flibe channels")
    mat_blanket_pre_tsq = get_material_by_name(Th_mats_pre_tsq, "doped flibe blanket")

    mat_channels_post_tsq = get_material_by_name(Th_mats_post_tsq, "doped flibe channels")
    mat_blanket_post_tsq = get_material_by_name(Th_mats_post_tsq, "doped flibe blanket")

    decay_heat_pre_tsq =  mat_channels_pre_tsq.get_decay_heat() + mat_blanket_pre_tsq.get_decay_heat()
    decay_heat_post_tsq =  mat_channels_post_tsq.get_decay_heat() + mat_blanket_post_tsq.get_decay_heat()

    if Th_fissile_masses[i, idx] < anp.sig_quantity:
        Th_decay_heats[i] = np.interp(Th_time_to_SQ[i], Th_time_steps[idx:idx+2], [decay_heat_pre_tsq, decay_heat_post_tsq])
    else:
        Th_decay_heats[i] = np.interp(Th_time_to_SQ[i], Th_time_steps[idx-1:idx+1], [decay_heat_pre_tsq, decay_heat_post_tsq])

    os.chdir('../../..')

print("Loaded decay heat data in "  + str(round(time.perf_counter() - init_time, 2)) + " seconds.")
# ====================================================
# Plotting
# ====================================================

masses = masses / 1e3 #convert from kg to metric tons

u_marker = 'o'
th_marker = 's'

u_color = 'tab:orange'
th_color = 'tab:purple'

dpi = 300

title_y = 1.05

fontdict = {"size":16}

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

ax.scatter(masses, U_time_to_SQ/24, label="U-238", marker=u_marker, color=u_color)
ax.scatter(masses, Th_time_to_SQ/24, label="Th-232", marker=th_marker, color=th_color)

ax.set_ylim(10, 200)

np.save("U_time_to_SQ_depletion", U_time_to_SQ)
np.save("Th_time_to_SQ_depletion", Th_time_to_SQ)

# Fit data to 1/x function:
def fit(x, A, B, C):
    return (A/x) - B*x + C

U_popt, U_pcov = curve_fit(fit, masses, U_time_to_SQ)
Th_popt, Th_pcov = curve_fit(fit, masses, Th_time_to_SQ)

fit_masses = np.linspace(1, masses[-1], num=100)
ax.plot(fit_masses, fit(fit_masses, *U_popt)/24, alpha=0.3, color=u_color)
ax.plot(fit_masses, fit(fit_masses, *Th_popt)/24, alpha=0.3, color=th_color)

ax.legend()

ax.set_xlim(0, masses[-1] + 2)
ax.set_ylim(10, np.max(Th_time_to_SQ/24) + 100)

ax.set_yscale("log")

ax.set_title("Time to Breed a Significant Quantity of Fissile Material", fontdict=fontdict, y=title_y)
ax.set_ylabel("Time (days)", fontdict=fontdict)
ax.set_xlabel("Mass of Fertile Material (metric tons)", fontdict=fontdict)

fig.savefig("time_to_sq.png", dpi=dpi)

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Fission Power

fig, ax = plt.subplots()
ax.spines["top"].set_color("None")
ax.spines["right"].set_color("None")

U_fission_power_at_SQ = np.empty((len(masses)))
Th_fission_power_at_SQ = np.empty((len(masses)))
for i, mass in enumerate(masses):
    """ Uranium """
    U_res = linregress(U_time_steps, U_fission_powers[i, :])
    U_fission_power_at_SQ[i] = U_res.intercept + U_res.slope*U_time_to_SQ[i]

    """ Thorium """
    Th_res = linregress(Th_time_steps, Th_fission_powers[i, ])
    Th_fission_power_at_SQ[i] = Th_res.intercept + Th_res.slope*Th_time_to_SQ[i]

size = 6

ax.scatter(masses, U_fission_powers[:, 0], marker=u_marker, c=u_color, s=size, label='U-238')
ax.scatter(masses, Th_fission_powers[:, 0], marker=th_marker, c=th_color, s=size, label='Th-232')

ax.scatter(masses, U_fission_power_at_SQ, marker=u_marker, c=u_color, s=size)
ax.scatter(masses, Th_fission_power_at_SQ, marker=th_marker, c=th_color, s=size)

ax.fill_between(masses, U_fission_powers[:, 0], U_fission_power_at_SQ, color=u_color, alpha=0.3)
ax.fill_between(masses, Th_fission_powers[:, 0], Th_fission_power_at_SQ, color=th_color, alpha=0.3)

text_offset = 5
ax.annotate("t = $t_{SQ}$", (masses[-1], U_fission_power_at_SQ[-1]), color=u_color, textcoords='offset points', xytext=(text_offset, 0))
ax.annotate("t = 0", (masses[-1], U_fission_powers[-1, 0]), color=u_color, textcoords='offset points', xytext=(text_offset, 0))

ax.annotate("t = $t_{SQ}$", (masses[-1], Th_fission_power_at_SQ[-1]), color=th_color, textcoords='offset points', xytext=(text_offset, 0))
ax.annotate("t = 0", (masses[-1], Th_fission_powers[-1, 0]), color=th_color, textcoords='offset points', xytext=(text_offset, 0))

ax.set_xlim(0, masses[-1] + 5)
ax.set_ylim(0, U_fission_power_at_SQ[-1] + 10)

ax.legend()

ax.set_title("Fission Power in Doped FLiBe Blanket", fontdict=fontdict, y=title_y)
ax.set_ylabel("Fission Power (MW)", fontdict=fontdict)
ax.set_xlabel("Fertile Mass (metric tons)", fontdict=fontdict)

fig.savefig("fission_power.png", dpi=dpi)

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Isotopic Purity

fig, ax = plt.subplots()
ax.spines["top"].set_color("None")
ax.spines["right"].set_color("None")

ax.scatter(masses, U_purities*100, label = "Pu-239", marker=u_marker, color=u_color)
ax.scatter(masses, Th_purities*100, label = "U-233", marker=th_marker, color=th_color)

ax.set_ylim(99.5, 100)

ax.legend()

ax.set_title("Isotopic Purity vs. Fertile Inventory", fontdict=fontdict, y=title_y)
ax.set_ylabel("Isotopic Purity (% fissile isotope)", fontdict=fontdict)
ax.set_xlabel("Fertile Mass (metric tons)", fontdict=fontdict)

fig.savefig("isotopic_purity.png", dpi=dpi)

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Flux Spectrum

fig, ax = plt.subplots()
ax.spines["top"].set_color("None")
ax.spines["right"].set_color("None")

for i, mass in enumerate(masses):
    ax.step(energies[1:], U_flux_spectra[i, -1, :], label=str(mass) + " kg")

ax.set_xlabel("Energy")
ax.set_ylabel("Flux (arb. units)")

ax.set_title("Average Neutron Flux Spectrum in Uranium Doped Blanket After 100 Days")

ax.set_ylim(0.01, 9)
ax.set_xlim(10, 1e8)

ax.set_xscale('log')
ax.set_yscale('log')

ax.legend()

fig.savefig("U_flux_spectra.png", dpi=300)

fig, ax = plt.subplots()
ax.spines["top"].set_color("None")
ax.spines["right"].set_color("None")

for i, mass in enumerate(masses):
    ax.step(energies[1:], Th_flux_spectra[i, -1, :], label=str(mass) + " kg")

ax.set_xlabel("Energy")
ax.set_ylabel("Flux (arb. units)")

ax.set_title("Average Neutron Flux Spectrum in Thorium Doped Blanket After 100 Days")

ax.set_ylim(0.01, 9)
ax.set_xlim(10, 1e8)

ax.set_xscale('log')
ax.set_yscale('log')

ax.legend()

fig.savefig("Th_flux_spectra.png", dpi=dpi)

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Flux Spectrum Evolution

# for i, mass in enumerate(masses):
#     fig, ax = plt.subplots()
#     ax.spines["top"].set_color("None")
#     ax.spines["right"].set_color("None")

#     energy_bin_centers = 0.5 * (energies[1:] + energies[:-1])

#     for j in range(0, num_steps):
#         difference_spectrum = (U_flux_spectra[i, j, :] - U_flux_spectra[i, 0, :])/U_flux_spectra[i, 0, :]

#         ax.step(energy_bin_centers, difference_spectrum, label=U_time_steps[j])

#     ax.set_yscale("log")
#     ax.set_xscale("log")
#     fig.savefig("U_spectrum_evolution_" + str(mass) + "_kg.png", dpi=300)

#     fig, ax = plt.subplots()
#     ax.spines["top"].set_color("None")
#     ax.spines["right"].set_color("None")

#     energy_bin_centers = 0.5 * (energies[1:] + energies[:-1])

#     for j in range(0, num_steps):
#         difference_spectrum = (Th_flux_spectra[i, j, :] - Th_flux_spectra[i, 0, :])/Th_flux_spectra[i, 0, :]

#         ax.step(energy_bin_centers, difference_spectrum, label=Th_time_steps[j])

#     ax.set_yscale("log")
#     ax.set_xscale("log")
#     fig.savefig("spectrum_evolution_" + str(mass) + "_kg.png", dpi=300)

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Decay Photon Spectra

# Uranium

fig, axs = plt.subplots(1, 2)

for ax in axs:
    ax.spines["top"].set_color("None")
    ax.spines["right"].set_color("None")

    ax.set_yscale("log")

    #ax.set_ylim(1e15, 1e19)
    #ax.set_xlim(0, 3)
    ax.set_xlabel("Photon Energy (MeV)")
    ax.set_ylabel("Bq/cm3/metric ton")

for i, mass in enumerate(masses):
    axs[0].scatter(U_decay_spectra_channels[i].x/1e6, U_decay_spectra_channels[i].p/mass, label = f"{mass} Metric Tons", s=1)
    axs[1].scatter(U_decay_spectra_blanket[i].x/1e6, U_decay_spectra_blanket[i].p/mass, label = f"{mass} Metric Tons", s=1)
    
axs[0].set_title("Gamma spectrum in a Uranium doped cooling channel")
axs[1].set_title("Gamma spectrum in a Uranium doped blanket")

#ax.legend()

fig.set_size_inches(10, 4)
fig.savefig("U_decay_spectra.png", dpi=dpi)

#Thorium:

fig, axs = plt.subplots(1, 2)

for ax in axs:
    ax.spines["top"].set_color("None")
    ax.spines["right"].set_color("None")

    ax.set_yscale("log")

    #ax.set_ylim(1e15, 1e19)
    #ax.set_xlim(0, 3)
    ax.set_xlabel("Photon Energy (MeV)")

for i in range(0, len(masses)):
    axs[0].step(Th_decay_spectra_channels[i].x/1e6, Th_decay_spectra_channels[i].p, label = str("Step " + str(int(i))))
    axs[1].step(Th_decay_spectra_blanket[i].x/1e6, Th_decay_spectra_blanket[i].p, label = str("Step " + str(int(i))))
    
axs[0].set_title("Gamma spectrum in a Thorium doped cooling channel")
axs[1].set_title("Gamma spectrum in a Thorium doped blanket")

fig.set_size_inches(10, 4)
fig.savefig("Th_decay_spectra.png", dpi=dpi)

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# TBR

fig, ax = plt.subplots()
ax.spines["top"].set_color("None")
ax.spines["right"].set_color("None")

ax.scatter(masses, U_TBR[:, 0, 0], marker=u_marker, color=u_color, label="U-238")
ax.scatter(masses, Th_TBR[:, 0, 0], marker=th_marker, color=th_color, label="Th-232")

#ax.plot(masses, U_TBR[:, 0, 0], color=u_color, alpha=0.3)
#ax.plot(masses, Th_TBR[:, 0, 0], color=th_color, alpha=0.3)

#ax.errorbar(masses, U_TBR[:, 1, 0], yerr=U_TBR[:, 1, 1], color="r", label="Uranium")
#ax.errorbar(masses, Th_TBR[:, 1, 0], yerr=Th_TBR[:, 1, 1], color="g", label="Thorium")

#ax.fill_between(masses, U_TBR[:, 0, 0], U_TBR[:, 1, 0], color='r', alpha=0.3)
#ax.fill_between(masses, Th_TBR[:, 0, 0], Th_TBR[:, 1, 0], color='g', alpha=0.3)

ax.set_ylim(1, 1.2)

ax.legend()

ax.set_title("TBR vs. Fertile Mass at $t=0$", fontdict=fontdict, y=title_y)
ax.set_ylabel("TBR", fontdict=fontdict)
ax.set_xlabel("Fertile Mass (metric tons)", fontdict=fontdict)

fig.savefig("fertile_tbr.png", dpi=dpi)

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Fissile Mass

for i, mass in enumerate(masses):
    fig, ax = plt.subplots()

    ax.spines["top"].set_color("None")
    ax.spines["right"].set_color("None")

    ax.plot(U_time_steps, U_fissile_masses[i], label="Pu239")
    ax.plot(Th_time_steps, Th_fissile_masses[i], label="U233")

    #U_fit = Polynomial.fit(U_time_steps, U_fissile_masses[i], 4)
    #Th_fit = Polynomial.fit(Th_time_steps, Th_fissile_masses[i], 4)

    #ax.plot(U_time_steps, U_fit(U_time_steps), alpha=0.5)
    #ax.plot(Th_time_steps, Th_fit(Th_time_steps), alpha=0.5)

    ax.hlines(anp.sig_quantity, Th_time_steps[0], Th_time_steps[-1], colors='tab:red', linestyles="dashed")

    ax.scatter(U_time_to_SQ[i]/24, anp.sig_quantity)
    ax.scatter(Th_time_to_SQ[i]/24, anp.sig_quantity)

    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Mass (kg)")
    ax.set_title("Fissile Mass vs. Time for a Fertile Mass of " + str(mass) + " metric tons")

    ax.legend()

    fig.savefig(str(mass) + "_metric_tons.png", dpi=dpi)

print("Completed Post Processing")

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Decay Heat

fig, ax = plt.subplots()

ax.spines["top"].set_color("None")
ax.spines["right"].set_color("None")

ax.scatter(masses, U_decay_heats/1e6, label="U-238", color=u_color, marker=u_marker)
ax.scatter(masses, Th_decay_heats/1e6, label="Th-238", color=th_color, marker=th_marker)

ax.set_ylim(0, 3.5)
ax.legend()

ax.set_title("Decay Heat vs. Fertile Mass at $t = t_{SQ}$", fontdict=fontdict)
ax.set_xlabel("Fertile Mass (Metric Tons)", fontdict=fontdict)
ax.set_ylabel("Decay Heat (MW)", fontdict=fontdict)

fig.savefig("decay_heat.png", dpi=dpi)

# ====================================================
# Data Storage
# ====================================================

U_data_dict ={"time_steps":U_time_steps,
              "time_to_sq":U_time_to_SQ,
              "fission_power_t_0":U_fission_powers[:, 0],
              "fission_power_t_sq":U_fission_power_at_SQ,
              "isotopic_purities":U_purities,
              "tbr_t0":U_TBR[:, 0, 0],
              "tbr_t_SQ":U_TBR[:, 1, 0],
              "flux_spectrum":U_flux_spectra,
              "fissile_mass":U_fissile_masses,
              "reaction_spectra":U_reaction_spectra,
              "decay_heat":U_decay_heats}

Th_data_dict ={"time_steps":Th_time_steps,
              "time_to_sq":Th_time_to_SQ,
              "fission_power_t_0":Th_fission_powers[:, 0],
              "fission_power_t_sq":Th_fission_power_at_SQ,
              "isotopic_purities":Th_purities,
              "tbr_t0":Th_TBR[:, 0, 0],
              "tbr_t_SQ":Th_TBR[:, 1, 0],
              "flux_spectrum":Th_flux_spectra,
              "fissile_mass":Th_fissile_masses,
              "reaction_spectra":Th_reaction_spectra,
              "decay_heat":Th_decay_heats}

os.chdir("../..")

try:
    os.chdir(base_dir + "/data")
except:
    os.mkdir(base_dir + "/data")
    os.chdir(base_dir + "/data")

with open("U_data_dict.pkl", 'wb') as file:
    pickle.dump(U_data_dict, file)

with open("Th_data_dict.pkl", 'wb') as file:
    pickle.dump(Th_data_dict, file)