import arc_nonproliferation as anp
import openmc
import openmc.deplete
import numpy as np
import os
import sys

""" Handle command line arguments """
if sys.argv[1] is not None:
    mass = int(sys.argv[1])
else:
    raise ValueError("Must specify a fertile mass in kg!")

if sys.argv[2] is not None:
    dopant = str(sys.argv[2])
else:
    raise ValueError("Must specify either U or Th as dopant type!")

# This function handles the simulation specific setup of each device object
def setup_device(device):
    """ Run settings """
    device.settings.photon_transport = False
    device.settings.particles = int(1e5)
    device.settings.batches = 10

    """ Cylindrical Mesh Tally """
    mesh = openmc.CylindricalMesh()
    mesh.r_grid = np.linspace(100, 700, num=50)
    mesh.z_grid = np.linspace(-300, 300, num=50)
    mesh.phi_grid = np.array([0, (2 * np.pi)/(18 * 2)])
    mesh_filter = openmc.MeshFilter(mesh)

    # """ Spectrum Mesh """
    # spec_mesh = openmc.CylindricalMesh()
    # spec_mesh.r_grid = np.linspace(100, 700, num=50)
    # spec_mesh.z_grid = np.linspace(-15, 15, num=1)
    # spec_mesh.phi_grid = np.array([0, (2 * np.pi)/(18 * 2)])
    # spec_mesh_filter = openmc.MeshFilter(mesh)

    """ Cell Filter """
    blanket_cell = device.get_cell(name='blanket')
    blanket_filter = openmc.CellFilter(blanket_cell)

    """ Energy Filter """
    energy_filter = openmc.EnergyFilter.from_group_structure("CCFE-709")

    device.add_tally('Mesh Tally', ['flux', '(n,Xt)', 'heating-local', 'absorption'], filters=[mesh_filter])

    """ FLiBe Tally """
    #flibe_filter = openmc.MaterialFilter(anp.get_material_by_name(device.materials, "doped_flibe"))
    device.add_tally('FLiBe Tally', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[])
    device.add_tally('Flux Tally', ['flux'], filters=[energy_filter, blanket_filter])
    device.add_tally('Li Tally', ['(n,Xt)'], filters=[], nuclides=['Li6', 'Li7'])
    device.add_tally("Feritle Tally", ['absorption', '(n,gamma)', 'fission'], nuclides=['U238, U235'])

    return device

# ==============================================================================
# Depletion Run
# ==============================================================================

mass = np.array([mass]) #kg of fertile material
np.savetxt('mass.txt', mass)

""" DEPLETION SETTINGS """
print("~~~~~~~~~~~~~~~~~~ FERTILE MASS: " + str(mass) + " kg ~~~~~~~~~~~~~~~~~~")

fusion_power = 500 #MW
num_steps = 10
time_steps = [15*24*60*60 / num_steps] * num_steps
source_rates = [fusion_power * anp.neutrons_per_MJ] * num_steps

chain_file = '/home/jlball/arc-nonproliferation/data/simple_fast_chain.xml'

""" Generate blankets doped to specified mass """

device = setup_device(anp.generate_device(dopant, mass[0]))

""" Run depletion calculation """

# This is necessary to ensure the tally.xml file gets written into the directory where the depletion calculation will run
device.build()
device.deplete(time_steps, 
    source_rates=source_rates, 
    method = 'epc_rk4',
    operator_kwargs={'chain_file':chain_file, 
                     'normalization_mode':'source-rate',
                     'dilute_initial':0, 
                     'reduce_chain':False,
                     'reduce_chain_level':3})

