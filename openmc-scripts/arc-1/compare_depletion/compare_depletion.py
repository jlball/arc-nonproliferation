import openmc
import openmc.deplete
import arc_nonproliferation.device as anp
import numpy as np
import os
import sys


###
"""
The purpose of this script is to run two depletion calculations with identical inputs, the only difference being once uses the 
independent operator and one uses the coupled operator.
"""
###


###
# DEPLETION SETTINGS
###
fusion_power = 500 #MW
num_steps = 10
time_steps = [100*24*60*60 / num_steps] * num_steps
source_rates = [fusion_power * anp.neutrons_per_MJ] * num_steps

chain_file = '/home/jlball/arc-nonproliferation/data/simple_chain_endfb71_pwr.xml'


# This function handles the simulation specific setup of each device object
def setup_device(device):
    """ Run settings """
    device.settings.photon_transport = False
    device.settings.particles = int(1e3)
    device.settings.batches = 10

    """ Cell Filter """
    blanket_cell = device.get_cell(name='blanket')
    blanket_filter = openmc.CellFilter(blanket_cell)

    """ Energy Filter """
    energy_filter = openmc.EnergyFilter.from_group_structure("CCFE-709")

    device.add_tally('FLiBe Tally', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[])
    device.add_tally('Flux Tally', ['flux'], filters=[energy_filter, blanket_filter])
    device.add_tally('Li Tally', ['(n,Xt)'], filters=[], nuclides=['Li6', 'Li7'])

    return device


# We now instantiate a model which has 20 metric tons of uranium dissolved in the blanket.
coupled_device = setup_device(anp.generate_device("U", 20e3))
independent_device = setup_device(anp.generate_device("U", 20e3))


### COUPLED DEPLETION ###
#Make a directory for and run the coupled depletion calculation

try:
    os.mkdir("COUPLED")
    os.chdir("COUPLED")
except:
    os.chdir("COUPLED")

coupled_device.build()

os.chdir('..')

coupled_device.deplete(time_steps, 
    source_rates=source_rates, 
    operator_kwargs={'chain_file':chain_file, 
                        'normalization_mode':'source-rate', 
                        'reduce_chain':False,
                        'reduce_chain_level':None}, 
    directory="COUPLED")

### INDEPENDENT DEPLETION ###
try:
    os.mkdir("INDEPENDENT")
    os.chdir("INDEPENDENT")
except:
    os.chdir("INDEPENDENT")

independent_device.build()

flux, micro_xs = openmc.deplete.get_microxs_and_flux(independent_device,
                                                        openmc.Materials.from_xml(),
                                                        chain_file=chain_file,
                                                        run_kwargs = {"threads":20,
                                                                    "particles":int(1e4)})

operator = openmc.deplete.IndependentOperator(openmc.Materials.from_xml(), 
                                                    flux,
                                                    micro_xs,
                                                    chain_file=chain_file, 
                                                    normalization_mode='source-rate', 
                                                    reduce_chain=False, 
                                                    reduce_chain_level=None, 
                                                    )

integrator = openmc.deplete.PredictorIntegrator(operator, 
                                                time_steps, 
                                                source_rates=source_rates,
                                                timestep_units='s')

integrator.integrate()