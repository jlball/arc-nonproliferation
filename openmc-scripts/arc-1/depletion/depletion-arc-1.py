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
    device.settings.particles = int(1e4)
    device.settings.batches = 10

    """ Cell Filter """
    blanket_cell = device.get_cell(name='blanket')
    blanket_filter = openmc.CellFilter(blanket_cell)

    """ Energy Filter """
    energy_filter = openmc.EnergyFilter.from_group_structure("CCFE-709")

    """ FLiBe Tally """
    #flibe_filter = openmc.MaterialFilter(anp.get_material_by_name(device.materials, "doped_flibe"))
    device.add_tally('FLiBe Tally', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[])
    device.add_tally('Flux Tally', ['flux'], filters=[energy_filter, blanket_filter])
    device.add_tally('Li Tally', ['(n,Xt)'], filters=[], nuclides=['Li6', 'Li7'])

    return device

# ==============================================================================
# Depletion Run
# ==============================================================================

mass = np.array([20e3]) #kg of fertile material
np.savetxt(base_dir + '/mass.txt', mass)

""" DEPLETION SETTINGS """
print("~~~~~~~~~~~~~~~~~~ FERTILE MASS: " + str(mass) + " kg ~~~~~~~~~~~~~~~~~~")

fusion_power = 500 #MW
num_steps = 50
time_steps = [100*24*60*60 / num_steps] * num_steps
source_rates = [fusion_power * anp.neutrons_per_MJ] * num_steps

chain_file = '/home/jlball/arc-nonproliferation/data/chain_endfb71_pwr.xml'

""" Generate blankets doped to specified mass """

U_device = setup_device(anp.generate_device("U", mass[0]))
Th_device = setup_device(anp.generate_device("Th", mass[0]))

""" Run depletion calculation """

# This is necessary to ensure the tally.xml file gets written into the directory where the depletion calculation will run
os.mkdir(base_dir + '/Uranium/'+ str(mass))
os.chdir(base_dir + '/Uranium/'+ str(mass))
U_device.build()

os.chdir('../../..')

U_device.deplete(time_steps, 
    source_rates=source_rates, 
    method = 'cecm',
    operator_kwargs={'chain_file':chain_file, 
                     'normalization_mode':'source-rate',
                     'reduce_chain':True,
                     'reduce_chain_level':5}, 
    directory=base_dir + '/Uranium/'+ str(mass))

os.mkdir(base_dir + '/Thorium/'+ str(mass))
os.chdir(base_dir + '/Thorium/'+ str(mass))
Th_device.build()
os.chdir('../../..')

Th_device.deplete(time_steps, 
    source_rates=source_rates, 
    method = 'cecm',
    operator_kwargs={'chain_file':chain_file, 
                     'normalization_mode':'source-rate',
                     'reduce_chain':True,
                     'reduce_chain_level':5}, 
    directory=base_dir + '/Thorium/' + str(mass))
