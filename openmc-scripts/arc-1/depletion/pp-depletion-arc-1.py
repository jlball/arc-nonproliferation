import openmc
from openmc.deplete import Results
import matplotlib.pyplot as plt
import numpy as np
import sys
from arc_nonproliferation.postprocess import *
from arc_nonproliferation.constants import *
from scipy.optimize import curve_fit
from scipy.stats import linregress
import time

if sys.argv[1] is not None:
    base_dir = './' + sys.argv[1]
else:
    raise ValueError("No base directory specified!")

fusion_power = 500 #MW
total_neutron_rate = fusion_power * neutrons_per_MJ

chain_file = '/home/jlball/arc-nonproliferation/data/simple_chain_endfb71_pwr.xml'
openmc.config['chain_file'] = chain_file

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Time to a Significant Quantity

""" Load masses and initialisze final output arrays """
mass = np.loadtxt(base_dir + '/mass.txt')
mass = np.array([float(mass)])
""" Extract time to 1 SQ for Uranium """
os.chdir(base_dir + "/Uranium/" + str(mass))

U_results = Results('depletion_results.h5')
U_time_to_SQ = extract_time_to_sq('U', U_results)

#While we're here, get the number of depletion steps:
time_steps = U_results.get_times()
num_steps = len(time_steps)

os.chdir("../../..")

""" Extract time to 1 SQ for Thorium """
os.chdir(base_dir + "/Thorium/" + str(mass))

Th_results = Results('depletion_results.h5')
Th_time_to_SQ = extract_time_to_sq('Th', Th_results)

os.chdir("../../..")

print("THORIUM TIME TO 1 SQ: " + str(Th_time_to_SQ/24) + " Days")
print("URANIUM TIME TO 1 SQ: " + str(U_time_to_SQ/24) + " Days")

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Fission Power

""" Uranium """
os.chdir(base_dir + "/Uranium/" + str(mass))

for step in range(0, num_steps):
    sp = openmc.StatePoint('openmc_simulation_n'+str(step)+'.h5')
    tally = sp.get_tally(name='FLiBe Tally')

    U_fission_power = anp.get_uvalue(tally, 'kappa-fission').n * total_neutron_rate * MJ_per_eV

os.chdir('../../..')

""" Thorium """
os.chdir(base_dir + "/Thorium/" + str(mass))

for step in range(0, num_steps):
    sp = openmc.StatePoint('openmc_simulation_n'+str(step)+'.h5')
    tally = sp.get_tally(name='FLiBe Tally')

    Th_fission_power = anp.get_uvalue(tally, 'kappa-fission').n * total_neutron_rate * MJ_per_eV

os.chdir('../../..')

print("Loaded fission power data...")

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Isotopic Purity

""" Uranium """
os.chdir(base_dir + "/Uranium/" + str(mass))

U_results = Results('depletion_results.h5')
U_purity = extract_isotopic_purity("U", U_results)
U_purity = U_purity[-1]

os.chdir('../../..')

""" Thorium """
os.chdir(base_dir + "/Thorium/" + str(mass))

Th_results = Results('depletion_results.h5')
Th_purity = extract_isotopic_purity("Th", Th_results)
Th_purity = Th_purity[-1]

os.chdir('../../..')

print("Loaded isotopic purity data...")

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Flux Spectrum in Blanket

U_flux_spectra = np.empty((num_steps, 709))
Th_flux_spectra = np.empty((num_steps, 709))

""" Uranium """
os.chdir(base_dir + "/Uranium/" + str(mass))

for step in range(0, num_steps):
    sp = openmc.StatePoint('openmc_simulation_n'+str(step)+'.h5')
    flux_tally = sp.get_tally(name='Flux Tally')
    flux_spectrum = flux_tally.get_reshaped_data()
    U_flux_spectra[step] = flux_spectrum.reshape((709,))

os.chdir('../../..')

""" Thorium """
os.chdir(base_dir + "/Thorium/" + str(mass))

for step in range(0, num_steps):
    sp = openmc.StatePoint('openmc_simulation_n'+str(step)+'.h5')
    flux_tally = sp.get_tally(name='Flux Tally')
    flux_spectrum = flux_tally.get_reshaped_data()
    Th_flux_spectra[step] = flux_spectrum.reshape((709,))

os.chdir('../../..')

energy_groups = openmc.mgxs.EnergyGroups(openmc.mgxs.GROUP_STRUCTURES['CCFE-709'])
energies = energy_groups.group_edges

print("Loaded flux spectrum data...")

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Absorption

U_absorption = np.empty((num_steps, 709, 3))
Th_absorption = np.empty((num_steps, 709, 3))

""" Uranium """
os.chdir(base_dir + "/Uranium/" + str(mass))

for step in range(0, num_steps):
    sp = openmc.StatePoint('openmc_simulation_n'+str(step)+'.h5')
    absorption_tally = sp.get_tally(name='Absorption Tally')
    absorption_spectrum = absorption_tally.get_reshaped_data()
    U_absorption[step] = absorption_spectrum.reshape((709,3))

os.chdir('../../..')

""" Thorium """
os.chdir(base_dir + "/Thorium/" + str(mass))

for step in range(0, num_steps):
    sp = openmc.StatePoint('openmc_simulation_n'+str(step)+'.h5')
    absorption_tally = sp.get_tally(name='Absorption Tally')
    absorption_spectrum = absorption_tally.get_reshaped_data()
    Th_absorption[step] = absorption_spectrum.reshape((709,3))

os.chdir('../../..')

print("Loaded absorption data...")


# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Decay Photon Spectrum

# U_decay_spectra = []
# Th_decay_spectra = []

# for i in range(0, num_steps):
#     """ Uranium """
#     os.chdir(base_dir + "/Uranium/" + str(mass))

#     U_results = Results('depletion_results.h5')
#     U_mats = U_results.export_to_materials(i)
#     flibe_mat = get_material_by_name(U_mats, "doped flibe blanket")
#     U_decay_spectra.append(flibe_mat.decay_photon_energy)

#     os.chdir('../../..')

#     """ Thorium """
#     os.chdir(base_dir + "/Thorium/" + str(mass))

#     Th_results = Results('depletion_results.h5')
#     Th_mats = Th_results.export_to_materials(i)
#     flibe_mat = get_material_by_name(Th_mats, "doped flibe blanket")
#     Th_decay_spectra.append(flibe_mat.decay_photon_energy)

#     os.chdir('../../..')

# print("Loaded decay photon spectrum data...")

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Fissile Mass
init_time = time.perf_counter()

""" Iterate through each mass simulated and get fissile mass at each time step"""
U_fissile_masses = np.empty((len(time_steps)))
Th_fissile_masses = np.empty((len(time_steps)))

os.chdir(base_dir + "/Uranium/" + str(mass))

U_results = Results('depletion_results.h5')
U_fissile_masses = get_masses_from_mats('Pu239', U_results)

os.chdir("../../..")

os.chdir(base_dir + "/Thorium/" + str(mass))

Th_results = Results('depletion_results.h5')
Th_fissile_masses = get_masses_from_mats('U233', Th_results)

os.chdir("../../..")

print("Loaded fissile mass data in "  + str(round(time.perf_counter() - init_time, 2)) + " seconds.")

# ====================================================
# Plotting
# ====================================================

# """ Uranium """
# os.chdir(base_dir + "/Uranium/" + str(mass))

# for step in range(0, num_steps):
#     sp = openmc.StatePoint('openmc_simulation_n'+str(step)+'.h5')
#     tally = sp.get_tally(name='Mesh Tally')

# fig, ax = anp.plot_RZ_quantity(mesh_tally, '(n,Xt)', volume_norm=False, title='Mesh TBR')

#Change into dedicated directory for figures or create figures directory
try:
    os.chdir(base_dir + "/figures")
except:
    os.mkdir(base_dir + "/figures")
    os.chdir(base_dir + "/figures")

fig, ax = plt.subplots()
ax.spines["top"].set_color("None")
ax.spines["right"].set_color("None")

for i, spectrum in enumerate(U_flux_spectra):
    ax.step(energies[1:], spectrum, label=str(i + 1))

ax.set_xlabel("Energy")
ax.set_ylabel("Flux (arb. units)")
ax.set_yscale('log')
ax.set_xscale('log')

ax.legend()
fig.savefig('U_flux_spectra.png')

fig, ax = plt.subplots()
ax.spines["top"].set_color("None")
ax.spines["right"].set_color("None")

ax.step(energies[1:], (U_flux_spectra[-1] - U_flux_spectra[0])*100/U_flux_spectra[0], label="Uranium")
ax.step(energies[1:], (Th_flux_spectra[-1] - Th_flux_spectra[0])*100/U_flux_spectra[0], label="Thorium")

ax.legend()

ax.set_xscale("log")

ax.set_xlabel("Energy")
ax.set_ylabel("Relative difference (percent)")
ax.set_ylim(0, 25)
fig.savefig('U_flux_spectra_difference.png')

# Decay Photon Spectrum

# fig, ax = plt.subplots()
# ax.spines["top"].set_color("None")
# ax.spines["right"].set_color("None")

# for i, dist in enumerate(U_decay_spectra):
#     ax.step(dist.x, dist.p, label = str("Step " + str(int(i))))

# ax.set_yscale("linear")

# ax.set_ylim(0, 1e17)

# ax.set_title("Gamma spectrum at $t_{SQ}$ in a Uranium doped blanket")
# ax.set_xlabel("Photon Energy (eV)")

# fig.savefig("U_decay_spectra.png")

# fig, ax = plt.subplots()
# ax.spines["top"].set_color("None")
# ax.spines["right"].set_color("None")

# for i, dist in enumerate(Th_decay_spectra):
#     ax.step(dist.x, dist.p, label = str("Step " + str(int(i))))

# ax.set_yscale("linear")

# ax.set_title("Gamma spectrum at $t_{SQ}$ in a Thorium doped blanket")
# ax.set_xlabel("Photon Energy (eV)")

# ax.set_ylim(0, 1e17)

# fig.savefig("Th_decay_spectra.png")

# Flux Spectrum Evolution

fig, ax = plt.subplots()
ax.spines["top"].set_color("None")
ax.spines["right"].set_color("None")

ax.set_yscale("log")

energy_bin_centers = 0.5 * (energies[1:] + energies[:-1])

for j in range(0, num_steps):
    difference_spectrum = (U_flux_spectra[j, :] - U_flux_spectra[0, :])/U_flux_spectra[0, :]

    ax.step(energy_bin_centers, difference_spectrum, label=time_steps[j])

fig.savefig("U_spectrum_evolution_" + str(mass) + "_kg.png", dpi=300)

fig, ax = plt.subplots()
ax.spines["top"].set_color("None")
ax.spines["right"].set_color("None")

ax.set_yscale('log')

energy_bin_centers = 0.5 * (energies[1:] + energies[:-1])

for j in range(0, num_steps):
    difference_spectrum = (Th_flux_spectra[j, :] - Th_flux_spectra[0, :])/Th_flux_spectra[0, :]

    ax.step(energy_bin_centers, difference_spectrum, label=time_steps[j])

fig.savefig("Th_spectrum_evolution_" + str(mass) + "_kg.png", dpi=300)


# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Fissile Mass

fig, ax = plt.subplots()

ax.plot(time_steps, U_fissile_masses, label="Pu239")
ax.plot(time_steps, Th_fissile_masses, label="U233")

ax.set_xlabel("Time (days)")
ax.set_ylabel("Mass (kg)")
ax.set_title("Fissile Mass vs. Time for a Fertile Mass of " + str(mass) + " metric tons")

ax.hlines(time_steps[0], time_steps[-1], [anp.sig_quantity])

ax.legend()

fig.savefig("fissile_mass.png", dpi=300)

# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
# Absorption

# First time step
fig, ax = plt.subplots()
ax.spines["top"].set_color("None")
ax.spines["right"].set_color("None")

ax.step(energy_bin_centers, U_absorption[0, :, -1], label="U238")
ax.step(energy_bin_centers, U_absorption[0, :, -1], label='Th232')

ax.set_xlabel("Energy")
ax.set_ylabel("Arb. Units")
ax.set_title("Absorption of neutrons by fertile isotopes at $t = 0$")

ax.set_xscale("log")
ax.set_yscale("log")

fig.savefig("Absorption_t0.png", dpi=300)

# Final timestep
fig, ax = plt.subplots()
ax.spines["top"].set_color("None")
ax.spines["right"].set_color("None")

ax.step(energy_bin_centers, U_absorption[-1 :, -1], label="U238")
ax.step(energy_bin_centers, U_absorption[-1 :, -1], label='Th232')

ax.set_xlabel("Energy")
ax.set_ylabel("Arb. Units")
ax.set_title("Absorption of neutrons by fertile isotopes at final time step")

ax.set_xscale("log")
ax.set_yscale("log")

fig.savefig("Absorption_tfinal.png", dpi=300)

