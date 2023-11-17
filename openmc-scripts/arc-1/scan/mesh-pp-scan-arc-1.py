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

""" Load masses and initialisze final output arrays """
masses = np.loadtxt(base_dir + '/masses.txt')
print(masses)
U_time_to_SQ = np.empty(len(masses))
Th_time_to_SQ = np.empty(len(masses))

# ====================================================
# TBR
# ====================================================

U_tbr = np.empty(len(masses))
Th_tbr = np.empty(len(masses))

for i, mass in enumerate(masses):

    """ Extract time to 1 SQ for Uranium """
    os.chdir(base_dir + "/Uranium/" + str(mass))
    sp = openmc.StatePoint('statepoint.10.h5')
    U_tally = sp.get_tally(name='Li Tally')

    U_r_mesh, U_z_mesh, U_abs_mesh = get_RZ_cyl_mesh_data(U_tally, "(n,gamma)", volume_norm=False)
    os.chdir("../../..")

    """ Extract time to 1 SQ for Thorium """
    os.chdir(base_dir + "/Thorium/" + str(mass))
    sp = openmc.StatePoint('statepoint.10.h5')
    Th_tally = sp.get_tally(name='Li Tally')

    Th_r_mesh, Th_z_mesh, Th_abs_mesh = get_RZ_cyl_mesh_data(Th_tally, "(n,gamma)", volume_norm=False)
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

fig, ax = plt.subplots()

cf = ax.contourf(U_r_mesh, U_z_mesh, U_abs_mesh)
ax.set_aspect(1)

fig.savefig("U_n_gamma_2D.png", dpi=300)