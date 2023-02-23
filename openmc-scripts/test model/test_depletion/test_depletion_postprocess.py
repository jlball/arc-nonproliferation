import openmc
from openmc.deplete import Results
import sys
import os
from arc_nonproliferation.postprocess import *
import arc_nonproliferation.constants as constants
import matplotlib.pyplot as plt

# TODO
# Clodgy fix, need to find a more elegant solution
openmc.config['chain_file'] = '/home/jlball/arc-nonproliferation/openmc-scripts/test_depletion/chain_endfb71_pwr.xml'

if sys.argv[1] is not None:
    base_dir = './' + sys.argv[1]
else:
    raise ValueError("No base directory specified!")

os.chdir(base_dir)

results = Results("depletion_results.h5")
time_steps = results.get_times('h')
num_steps = len(time_steps)
 
fusion_power = 500 #MW
total_source_rate = fusion_power * constants.neutrons_per_MJ # neutrons/s

# ==============================================================================
# Fissile Inventory
# ==============================================================================
U_fissile_masses = get_masses_from_mats("U", results)

# ==============================================================================
# Fission Power
# ==============================================================================
fission_powers = np.empty(num_steps)
fission_powers_errs = np.empty(num_steps)
for step in range(0, num_steps):
    sp = openmc.StatePoint('openmc_simulation_n' + str(step) + '.h5')

    tally = sp.get_tally(name='FLiBe Tally')
    fission_powers[step] = get_uvalue(tally, 'kappa-fission').n * total_source_rate * constants.MJ_per_eV
    fission_powers_errs[step] = get_uvalue(tally, 'kappa-fission').s * total_source_rate * constants.MJ_per_eV

# ==============================================================================
# Decay Heat
# ==============================================================================
decay_heats = extract_decay_heat(results)

# ==============================================================================
# Gamma Spectrum
# ==============================================================================
# Currently only displays gamma spectrum at LAST STEP
materials = results.export_to_materials(num_steps - 1)
doped_flibe = get_material_by_name(materials, 'doped flibe') 
gamma_energies = doped_flibe.decay_photon_energy
#gamma_intensities = doped_flibe.decay_

# ==============================================================================
# Plotting
# ==============================================================================

#Change into dedicated directory for figures or create figures directory
try:
    os.chdir("figures")
except:
    os.mkdir("figures")
    os.chdir("figures")

""" Fissile Mass """
fig, ax = plt.subplots()
ax.scatter(time_steps, U_fissile_masses, label="Pu-239")
ax.legend()

ax.set_title("Mass of Pu-239 vs. Time")
ax.set_xlabel("Time (hours)")
ax.set_ylabel("Mass of Fissile Material (kg)")

fig.savefig("fissile_mass.png")

""" Fission Power """
fig, ax = plt.subplots()
ax.errorbar(time_steps, fission_powers, yerr=fission_powers_errs)

ax.set_title("Total Fission Power vs. Time")
ax.set_xlabel("Time (hours)")
ax.set_ylabel("Fission Power (MW)")

fig.savefig("fission_power.png")

""" Gamma Spectrum """
fig, ax = plt.subplots()
#ax.bar(gamma_energies, gamma_intensities)

ax.set_title("Gamma Spectrum in Blanket at Last Timestep")
ax.set_xlabel("Photon Energy (eV)")
ax.set_ylabel("Intensity (arb. units)")

fig.savefig("gamma_spectrum.png")

""" Decay Heat """
fig, ax = plt.subplots()
ax.scatter(time_steps, decay_heats)

ax.set_title("Decay Heat in Blanket vs. Time")
ax.set_xlabel("Time (hours)")
ax.set_ylabel("Power (MW)")
ax.set_yscale("log")

fig.savefig("decay_heat.png")
