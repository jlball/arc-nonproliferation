import openmc
from dose_rate_model import generate_dose_rate_model, dose_tally_volume
from openmc.deplete import Results
import sys
import os
import numpy as np
import pickle
from arc_nonproliferation.postprocess import *
from arc_nonproliferation.device import *
from arc_nonproliferation.constants import chain_file
from arc_nonproliferation.materials import make_doped_flibe, air
from decay_calc import decay_time_steps

dose_rate_folder_name = "dose_rate_calc"
cooldown_folder_name = "dose_rate_cooldown"

openmc.config['chain_file'] = chain_file

dopants = ["U", "Th"]

# Get volume of blanket and channel cells
vol_calc_load = openmc.VolumeCalculation.from_hdf5('/home/jlball/arc-nonproliferation/data/arc-1_volumes.h5')
blanket_volume = vol_calc_load.volumes[8].n
channels_volume = vol_calc_load.volumes[5].n

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

    dose_rates = np.zeros((len(masses), len(decay_time_steps)+1))
    for i, mass in enumerate(masses):
        os.chdir(f"{mass}")
        os.chdir(f"{cooldown_folder_name}")

        results = Results("depletion_results.h5")

        time_steps = results.get_times()

        for j, result in enumerate(results):

            print(f"TIME STEP {j}: {time_steps[j]}")
            mats = results.export_to_materials(j)

            # Linearly interpolate material compositions to t_SQ
            blanket_mat = result.get_material(str(get_material_by_name(mats, "doped flibe blanket").id))
            channel_mat = result.get_material(str(get_material_by_name(mats, "doped flibe channels").id))

            blanket_mat.volume = blanket_volume
            channel_mat.volume = channels_volume

            blanket_mat = cutoff_nuclides(blanket_mat, 1e-200)
            channel_mat = cutoff_nuclides(channel_mat, 1e-200)

            # Create or enter subdirectory for decay calc
            try:
                os.chdir(f"{dose_rate_folder_name}_{j}")
            except:
                os.mkdir(f"{dose_rate_folder_name}_{j}")
                os.chdir(f"{dose_rate_folder_name}_{j}")

            # Perform decay only depletion calc
            model, dose_tally_volume = generate_dose_rate_model(blanket_mat, channel_mat, make_doped_flibe(dopant, mass, volume=blanket_volume+channels_volume))
            sp_path = model.run()

            # Postprocess to get dose rate in Sv/hr
            sp = openmc.StatePoint(sp_path)

            tally = sp.get_tally(name="dose_tally")
            dose_rate = (tally.get_values()[0][0][0])/dose_tally_volume # Normalize by the tally volume
            dose_rate = dose_rate * 1E-12 # convert from pSV/hr to Sv/hr
            dose_rates[i, j] = dose_rate
            os.chdir("..")

        os.chdir("../..")

    data_dict["dose_rate_cooldown"] = dose_rates
    data_dict["dose_rate_cooldown_times"] = time_steps

    os.chdir("../..")

    with open(f'{base_dir}/data/{dopant}_data_dict.pkl', 'wb') as file:
        pickle.dump(data_dict, file)

    
    