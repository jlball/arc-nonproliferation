import arc_nonproliferation as anp
import openmc
import openmc.deplete
import numpy as np
import os
import sys
import pickle

"""
This script runs an independent depletion calculation for the arc-1 model for a 
set of blanket materials with varying fertile masses present. Used to asses the
time needed to breed a significant quantity of fissile material. The output of 
this scrit can be automatically analysed using the corresponding script:

pp-independent-depletion-arc-1.py

Usage: Specify a set of fertile masses you wish to simulate, and a set of
desired run parameters for each calculation of fluxes and micro_xs in the 
setup_device function. Fusion power and timesteps can be adjusted in the 
depletion run section. To run the script, launch it from the command line
with a single additional argument specifying the name of the directory in 
which the output will be stored:

python3 independent-depletion-arc-1.py <output_directory_name>

"""

# ==============================================================================
# Setup
# ==============================================================================

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
    device.settings.particles = int(1e5)
    device.settings.batches = 10
    device.settings.survival_biasing = True

    return device

# ==============================================================================
# Depletion Run
# ==============================================================================

masses = np.array([5e3, 7e3, 10e3, 15e3, 30e3, 50e3]) #kg of fertile material
np.savetxt(base_dir + '/masses.txt', masses) # Store masses for later use in post processing


""" DEPLETION SETTINGS """
num_steps = 50
time_steps = [365 / num_steps] * num_steps
source_rates = [anp.Device.neutron_source_rate] * num_steps

openmc.config['chain_file'] = anp.constants.chain_file

for mass in masses:
    """ Generate blankets doped to specified mass """

    U_device = setup_device(anp.generate_device("U", mass))
    Th_device = setup_device(anp.generate_device("Th", mass))

    """ Run Uranium depletion calculation """
    print("~~~~~~~~~~~~~~~~~~ FERTILE MASS: " + str(mass) + " kg" + " DOPANT: Uranium-238" + " ~~~~~~~~~~~~~~~~~~")

    # This is necessary to ensure the tally.xml file gets written into the directory where the depletion calculation will run
    os.mkdir(base_dir + '/Uranium/'+ str(mass))
    os.chdir(base_dir + '/Uranium/'+ str(mass))
    U_device.build()

    # Generate fluxes and microxs from model
    U_flux, U_micro_xs = openmc.deplete.get_microxs_and_flux(U_device,
                                                            [U_device.channel, U_device.blanket],
                                                            run_kwargs = {"threads":20,
                                                                        "particles":int(1e3)})

    # Save flux and microxs to disk so transport calculation does not need to be rerun in future
    flux_file = open('U_flux', 'ab')
    pickle.dump(U_flux, flux_file)

    microxs_file = open("U_microxs", 'ab')
    pickle.dump(U_micro_xs, microxs_file)

    # Run depletion calculation
    U_operator = openmc.deplete.IndependentOperator(openmc.Materials([U_device.channel.fill, U_device.blanket.fill]), 
                                                    U_flux,
                                                    U_micro_xs,
                                                    normalization_mode='source-rate', 
                                                    reduce_chain=True, 
                                                    reduce_chain_level=5, 
                                                    )

    U_integrator = openmc.deplete.PredictorIntegrator(U_operator, 
                                                    time_steps, 
                                                    source_rates=source_rates,
                                                    timestep_units='d')

    U_integrator.integrate()

    os.chdir('../../..')

    """ Run Thorium depletion calculation """
    print("~~~~~~~~~~~~~~~~~~ FERTILE MASS: " + str(mass) + " kg" + " DOPANT: Thorium-232" + " ~~~~~~~~~~~~~~~~~~")

     # This is necessary to ensure the tally.xml file gets written into the directory where the depletion calculation will run
    os.mkdir(base_dir + '/Thorium/'+ str(mass))
    os.chdir(base_dir + '/Thorium/'+ str(mass))
    Th_device.build()

    # Generate fluxes and microxs from model
    Th_flux, Th_micro_xs = openmc.deplete.get_microxs_and_flux(Th_device,
                                                                    [Th_device.channel, Th_device.blanket],
                                                                    run_kwargs = {"threads":20,
                                                                                    "particles":int(1e3)})

    # Save flux and microxs to disk so transport calculation does not need to be rerun in future
    flux_file = open('U_flux', 'ab')
    pickle.dump(U_flux, flux_file)

    microxs_file = open("U_microxs", 'ab')
    pickle.dump(U_micro_xs, microxs_file)

    # Run depletion calculation
    Th_operator = openmc.deplete.IndependentOperator(openmc.Materials([Th_device.channel.fill, Th_device.blanket.fill]), 
                                                    Th_flux,
                                                    Th_micro_xs,
                                                    normalization_mode='source-rate', 
                                                    reduce_chain=True, 
                                                    reduce_chain_level=5, 
                                                    )

    Th_integrator = openmc.deplete.PredictorIntegrator(Th_operator, 
                                                    time_steps, 
                                                    source_rates=source_rates,
                                                    timestep_units='d')

    Th_integrator.integrate()

    os.chdir('../../..')

