import openmc
from openmc.deplete import Results, IndependentOperator, PredictorIntegrator
import sys
import os
import numpy as np
import pickle
from arc_nonproliferation.postprocess import *
from arc_nonproliferation.device import *
from arc_nonproliferation.constants import chain_file

cooldow_folder_name = "dose_rate_cooldown"

decay_time_steps = np.logspace(-1, 4, num=8) #days
source_rates = np.zeros(len(decay_time_steps))

openmc.config['chain_file'] = chain_file

dopants = ["U", "Th"]

if sys.argv[1] is not None:
    base_dir = './' + sys.argv[1]
else:
    raise ValueError("No base directory specified!")

masses = np.loadtxt(base_dir + "/masses.txt")

for dopant in dopants:
    with open(f'{base_dir}/data/{dopant}_data_dict.pkl', 'rb') as file:
        data_dict = pickle.load(file)

        time_to_sq = data_dict["time_to_sq"]/24

    if dopant == "U":
        os.chdir(f"{base_dir}/Uranium")
    elif dopant == "Th":
        os.chdir(f"{base_dir}/Thorium")
    else:
        raise ValueError("Invalid dopant type")

    for i, mass in enumerate(masses):
        os.chdir(f"{mass}")

        results = Results("depletion_results.h5")

        # Linearly interpolate material compositions to t_SQ
        blanket_mat = lin_interp_material(results, "doped flibe blanket", time_to_sq[i])
        channel_mat = lin_interp_material(results, "doped flibe channels", time_to_sq[i])

        # Create or enter subdirectory for decay calc
        try:
            os.chdir(cooldow_folder_name)
        except:
            os.mkdir(cooldow_folder_name)
            os.chdir(cooldow_folder_name)

        # Perform decay only depletion calc

        device = generate_device("U", 0, Li6_enrichment=7.5)

        # Set new material volumes
        channel_mat.volume = device.channel.fill.volume
        blanket_mat.volume = device.channel.fill.volume

        channel_mat.name =  device.channel.fill.name
        blanket_mat.name = device.blanket.fill.name


        device.channel.fill = channel_mat
        device.blanket.fill = blanket_mat

        device.build()

        operator = IndependentOperator(openmc.Materials([device.channel.fill, device.blanket.fill]),
                                        fluxes=np.zeros(2),
                                        micros=[openmc.deplete.MicroXS(np.zeros((0,0)), [], []), openmc.deplete.MicroXS(np.zeros((0,0)), [], [])],
                                        normalization_mode="source-rate")

        integrator = PredictorIntegrator(operator, 
                                            decay_time_steps, 
                                            source_rates=source_rates,
                                            timestep_units='d')
        
        integrator.integrate()

        os.chdir("../..")
    
    os.chdir("../..")
