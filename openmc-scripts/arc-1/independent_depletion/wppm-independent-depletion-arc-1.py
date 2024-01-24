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
    device.settings.survival_biasing = True

    return device

# ==============================================================================
# Depletion Run
# ==============================================================================

wppms = np.array([50, 100, 150, 200]) #kg of fertile material
np.savetxt(base_dir + '/masses.txt', wppms) # Store masses for later use in post processing

""" DEPLETION SETTINGS """
num_steps = 50
time_steps = [30 *365 / num_steps] * num_steps
source_rates = [anp.Device.neutron_source_rate] * num_steps

openmc.config['chain_file'] = anp.constants.chain_file

for wppm in wppms:
    """ Generate blankets doped to specified mass """

    U_device = setup_device(anp.generate_device("U", wppm, Li6_enrichment = Li6_enrichment, dopant_mass_units="wppm"))

    """ Run Uranium depletion calculation """
    print(f"~~~~~~~~~~~~~~~~~~ FERTILE MASS: {wppm} wppm DOPANT: Uranium-238 ~~~~~~~~~~~~~~~~~~")

    # This is necessary to ensure the tally.xml file gets written into the directory where the depletion calculation will run
    os.mkdir(base_dir + '/Uranium/'+ str(wppm))
    os.chdir(base_dir + '/Uranium/'+ str(wppm))
    U_device.build()

    # Generate fluxes and microxs from model
    U_flux, U_micro_xs = openmc.deplete.get_microxs_and_flux(U_device,
                                                            [U_device.channel, U_device.blanket])

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

