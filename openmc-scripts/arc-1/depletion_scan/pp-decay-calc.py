import openmc
from openmc.deplete import Results
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import numpy as np
import sys
from arc_nonproliferation.postprocess import *
from arc_nonproliferation.constants import *
from progressbar import progressbar
import pickle
import multiprocessing as mp

cooldow_folder_name = "dose_rate_cooldown"
num_timesteps = 51

openmc.config['chain_file'] = chain_file

dopants = ["U", "Th"]

use_stored_data = False

if sys.argv[1] is not None:
    base_dir = './' + sys.argv[1]
else:
    raise ValueError("No base directory specified!")

masses = np.loadtxt(base_dir + "/masses.txt")

dose_rate_dict = {}

if __name__ == "__main__":
    if not use_stored_data:
        os.chdir(f"{base_dir}")

        U_dose_rates, Th_dose_rates = np.zeros((len(masses), num_timesteps)), np.zeros((len(masses), num_timesteps))

        pool = mp.Pool()

        U_blanket_args = []
        U_channel_args = []

        for i, mass in enumerate(masses):
            print(f"~~~~~~ MASS: {mass} ~~~~~~~")

            # Get blanket and channel matertials at each timestep for Uranium
            os.chdir(f"Uranium/{mass}/{cooldow_folder_name}")
            U_results = Results("depletion_results.h5")
            time_steps = U_results.get_times()

            U_blanket_mats = []
            U_channel_mats = []
            for j, time in enumerate(time_steps):
                U_blanket_mats.append(get_material_by_name(U_results.export_to_materials(j), "doped flibe blanket"))
                U_channel_mats.append(get_material_by_name(U_results.export_to_materials(j), "doped flibe channels"))
            os.chdir("../../..")

            # Get blanket and channel matertials at each timestep for Thorium
            os.chdir(f"Thorium/{mass}/{cooldow_folder_name}")
            Th_results = Results("depletion_results.h5")

            Th_blanket_mats = []
            Th_channel_mats = []
            for j, time in enumerate(time_steps):
                Th_blanket_mats.append(get_material_by_name(Th_results.export_to_materials(j), "doped flibe blanket"))
                Th_channel_mats.append(get_material_by_name(Th_results.export_to_materials(j), "doped flibe channels"))
            os.chdir("../../..")

            # Dispatch computation of contact dose rates to process pool for speed
            U_blanket_result = pool.map_async(extract_contact_dose_rate, U_blanket_mats)
            U_channel_result = pool.map_async(extract_contact_dose_rate, U_channel_mats)

            Th_blanket_result = pool.map_async(extract_contact_dose_rate, Th_blanket_mats)
            Th_channel_result = pool.map_async(extract_contact_dose_rate, Th_channel_mats)

            # Sum results element-wise and store in 2D array
            U_dose_rates[i] = np.asarray(U_blanket_result.get()) + np.asarray(U_channel_result.get())
            Th_dose_rates[i] = np.asarray(Th_blanket_result.get()) + np.asarray(Th_channel_result.get())

            if mass == masses[2]:
                with open("U_blanket_mats.pkl", 'wb') as file:
                    pickle.dump(U_blanket_mats, file)

        dose_rate_dict["U"] = U_dose_rates
        dose_rate_dict["Th"] = Th_dose_rates
        dose_rate_dict["time_steps"] = time_steps

        # Store dose rate data in a pickle file:
        try:
            os.chdir("data")
        except:
            os.mkdir("data")
            os.chdir("data")

        with open("dose_rate_dict.pkl", 'wb') as file:
            pickle.dump(dose_rate_dict, file)

    else:
        with open(f"{base_dir}/data/dose_rate_dict.pkl", 'rb') as file:
            dose_rate_dict = pickle.load(file)

    try:
        os.chdir(f"{base_dir}/figures")
    except:
        os.mkdir(f"{base_dir}/figures")
        os.chdir(f"{base_dir}/figures")

    for dopant in dopants:
        dose_rates = dose_rate_dict[dopant]
        time_steps = dose_rate_dict["time_steps"]

        if dopant == 'U':
            plt_color ='r'
            plt_cm = cm.Oranges
        else:
            plt_color = 'g'
            plt_cm = cm.Purples

        norm = colors.Normalize(vmin=-0.5*masses[-1], vmax=1.1*masses[-1])

        fig, ax = plt.subplots()

        for i, dose_rate in enumerate(dose_rates):
            ax.plot(time_steps, dose_rates[i], color=plt_cm(norm(masses[i])), label=int(masses[i]/1e3))

        ax.set_xlabel("Time (days)")
        ax.set_ylabel("Contact Dose Rate (Sv/hr)")

        if dopant == "U":
            ax.set_title("Decay of contact dose rate in Uranium fueled FLiBe")

        if dopant == "Th":
            ax.set_title("Decay of contact dose rate in Thorium fueled FLiBe")

        ax.set_xscale("log")
        ax.legend()

        fig.savefig(f"{dopant}_contact_dose_rate_decay.png", dpi=300)

    os.chdir("../..")



