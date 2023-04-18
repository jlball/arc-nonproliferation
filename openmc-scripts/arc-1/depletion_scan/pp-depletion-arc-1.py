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

# ====================================================
# Fission Power
# ====================================================

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

# ====================================================
# Mesh Tally
# ====================================================



# ====================================================
# Isotopic Purity
# ====================================================
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

# ====================================================
# Flux Spectrum in Blanket
# ====================================================

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

