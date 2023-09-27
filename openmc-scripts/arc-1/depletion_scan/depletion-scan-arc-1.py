import arc_nonproliferation.device as anp
import arc_nonproliferation.constants as anp_constants
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

""" allows for enrichment to specified, defaults to natual """
if sys.argv[2] is not None:
    Li6_enrichment = float(sys.argv[2])
else:
    Li6_enrichment = 7.5

# This function handles the simulation specific setup of each device object
def setup_device(device):
    """ Run settings """
    device.settings.photon_transport = False
    device.settings.particles = int(1e4)
    device.settings.batches = 100
    device.survival_biasing = True

    """ Cell Filter """
    blanket_cell = device.get_cell(name='blanket')
    blanket_filter = openmc.CellFilter(blanket_cell)

    channel_cell = device.get_cell(name='channel')
    channel_filter = openmc.CellFilter(channel_cell)

    """ Energy Filter """
    energy_filter = openmc.EnergyFilter.from_group_structure("CCFE-709")

    """ FLiBe Tally """
    device.add_tally('FLiBe Tally', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[])
    device.add_tally('Flux Tally', ['flux'], filters=[energy_filter, blanket_filter, channel_filter])
    device.add_tally('Li Tally', ['(n,Xt)'], filters=[energy_filter, blanket_filter, channel_filter], nuclides=['Li6', 'Li7'])

    """ Fertile/Fissile Tally """
    if device.dopant == "U":
        fertile_nuclide = "U238"
        fissile_nuclide = "Pu239"
    elif device.dopant == "Th":
        fertile_nuclide = "Th232"
        fissile_nuclide = "U233"
    else:
        raise ValueError("Invalid Dopant Type!")

    device.add_tally("Fertile Tally", ['fission','(n,gamma)'], filters=[energy_filter, blanket_filter, channel_filter], nuclides=[fertile_nuclide])
    device.add_tally('Fissile Tally', ['fission', '(n,gamma)'], filters=[energy_filter, blanket_filter, channel_filter], nuclides=[fissile_nuclide])

    return device

# ==============================================================================
# Depletion Scan
# ==============================================================================

openmc.config['chain_file'] = anp_constants.chain_file

masses = np.array([5e3, 10e3, 20e3, 30e3, 40e3, 50e3])

np.savetxt(base_dir + '/masses.txt', masses)

for mass in masses:
    """ DEPLETION SETTINGS """
    print("~~~~~~~~~~~~~~~~~~ FERTILE MASS: " + str(mass) + " kg ~~~~~~~~~~~~~~~~~~")

    fusion_power = 500 #MW
    decay_num_steps = 9 # Number of timesteps to take within 3 half lives of breeding decay
    linear_num_steps = 10 # Number of timesteps to take beyond 3 half lives of breeding decay

    linear_time_steps = np.array([100*24*60*60 / linear_num_steps] * linear_num_steps)

    #Generate first set of timesteps based on decay of logest life isotope in breeding decay chain
    U_time_steps = np.linspace(0, 3/openmc.data.decay_constant("Np239"), num=decay_num_steps+1)

    U_time_steps = U_time_steps[1:] - U_time_steps[:-1]
    U_time_steps = np.append(U_time_steps, linear_time_steps)

    Th_time_steps = np.linspace(0, 3/openmc.data.decay_constant('Pa233'), num=decay_num_steps + 1)

    Th_time_steps = Th_time_steps[1:] - Th_time_steps[:-1]
    Th_time_steps = np.append(Th_time_steps, linear_time_steps)

    # Setup constant array of source rates
    source_rates = [fusion_power * anp.neutrons_per_MJ] * (decay_num_steps + linear_num_steps)

    """ Generate blankets doped to specified mass """
    
    U_device = setup_device(anp.generate_device("U", mass))
    Th_device = setup_device(anp.generate_device("Th", mass))

    """ Run depletion calculation """

    # This is necessary to ensure the tally.xml file gets written into the directory where the depletion calculation will run
    os.mkdir(base_dir + '/Uranium/'+ str(mass))
    os.chdir(base_dir + '/Uranium/'+ str(mass))
    U_device.build()

    os.chdir('../../..')

    U_device.deplete(U_time_steps, 
        source_rates=source_rates, 
        operator_kwargs={'normalization_mode':'source-rate', 
                         'reduce_chain':True,
                         'reduce_chain_level':5}, 
        directory=base_dir + '/Uranium/'+ str(mass),
        timestep_units='s',
        method='cecm')

    os.mkdir(base_dir + '/Thorium/'+ str(mass))
    os.chdir(base_dir + '/Thorium/'+ str(mass))
    Th_device.build()
    os.chdir('../../..')

    Th_device.deplete(Th_time_steps, 
        source_rates=source_rates, 
        operator_kwargs={'normalization_mode':'source-rate',
                         'reduce_chain':True,
                         'reduce_chain_level':5}, 
        directory=base_dir + '/Thorium/' + str(mass),
        timestep_units='s',
        method='cecm')
