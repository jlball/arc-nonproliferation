import arc_nonproliferation as anp
import openmc
import openmc.deplete
import numpy as np
import os
import sys

""" Handle command line arguments and setup directory structure """
if sys.argv[1] is not None:
    base_dir = str(sys.argv[1])
    os.mkdir(base_dir)
    print("OpenMC output saved to:", base_dir)

    os.mkdir(base_dir + '/Uranium')
    os.mkdir(base_dir + '/Thorium')

# This function handles the simulation specific setup of each device object
def setup_device(device):
    """ Run settings """
    device.settings.photon_transport = False
    device.settings.particles = int(1e3)
    device.settings.batches = 10

    """ Cylindrical Mesh Tally """
    mesh = openmc.CylindricalMesh()
    mesh.r_grid = np.linspace(100, 700, num=50)
    mesh.z_grid = np.linspace(-300, 300, num=50)
    mesh.phi_grid = np.array([0, (2 * np.pi)/(18 * 2)])
    mesh_filter = openmc.MeshFilter(mesh)

    device.add_tally('Mesh Tally', ['flux', '(n,Xt)', 'heating-local', 'absorption'], filters=[mesh_filter])

    """ Cell Filter """
    blanket_cell = device.get_cell(name='blanket')
    blanket_filter = openmc.CellFilter(blanket_cell)

    """ Energy Filter """
    energy_filter = openmc.EnergyFilter.from_group_structure("CCFE-709")

    """ FLiBe Tally """
    #flibe_filter = openmc.MaterialFilter(anp.get_material_by_name(device.materials, "doped_flibe"))
    
    # if device.dopant == "U":
    #     device.add_tally('Fissile Fission Tally', ['fission'], filters=[], nuclides=["Pu239"])

    # elif device.dopant == "Th":
    #     device.add_tally('Fissile Fission Tally', ['fission'], filters=[], nuclides=["U233"])

    device.add_tally('FLiBe Tally', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[])
    device.add_tally('Flux Tally', ['flux'], filters=[energy_filter, blanket_filter])
    device.add_tally('Li Tally', ['(n,Xt)'], filters=[], nuclides=['Li6', 'Li7'])

    return device

# ==============================================================================
# Depletion Scan
# ==============================================================================

masses = np.array([5e3, 10e3])

np.savetxt(base_dir + '/masses.txt', masses)

for mass in masses:
    """ DEPLETION SETTINGS """
    print("~~~~~~~~~~~~~~~~~~ FERTILE MASS: " + str(mass) + " kg ~~~~~~~~~~~~~~~~~~")

    fusion_power = 500 #MW
    num_steps = 3
    time_steps = [365*24*60*60 / num_steps] * num_steps
    source_rates = [fusion_power * anp.neutrons_per_MJ] * num_steps

    chain_file = '/home/jlball/arc-nonproliferation/data/simple_chain_endfb71_pwr.xml'

    """ Generate blankets doped to specified mass """
    
    U_device = setup_device(anp.generate_device("U", mass))
    Th_device = setup_device(anp.generate_device("Th", mass))

    """ Run depletion calculation """

    # This is necessary to ensure the tally.xml file gets written into the directory where the depletion calculation will run
    os.mkdir(base_dir + '/Uranium/'+ str(mass))
    os.chdir(base_dir + '/Uranium/'+ str(mass))
    U_device.build()

    os.chdir('../../..')

    U_device.deplete(time_steps, 
        source_rates=source_rates, 
        operator_kwargs={'chain_file':chain_file, 
                         'normalization_mode':'source-rate', 
                         'dilute_initial':0, 
                         'reduce_chain':False,
                         'reduce_chain_level':3}, 
        directory=base_dir + '/Uranium/'+ str(mass))

    os.mkdir(base_dir + '/Thorium/'+ str(mass))
    os.chdir(base_dir + '/Thorium/'+ str(mass))
    Th_device.build()
    os.chdir('../../..')

    Th_device.deplete(time_steps, 
        source_rates=source_rates, 
        operator_kwargs={'chain_file':chain_file, 
                         'normalization_mode':'source-rate',
                         'dilute_initial':0, 
                         'reduce_chain':False,
                         'reduce_chain_level':3}, 
        directory=base_dir + '/Thorium/' + str(mass))
