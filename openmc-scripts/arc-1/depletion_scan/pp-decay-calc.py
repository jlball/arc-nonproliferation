import openmc
from openmc.deplete import Results
import matplotlib.pyplot as plt
import numpy as np
import sys
from arc_nonproliferation.postprocess import *
from arc_nonproliferation.constants import *
from progressbar import progressbar
import pickle
from multiprocessing import Process

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


def parallel_extract_dose_rate(all_mats, material_name, time_steps, output_variable):
    dose_rates = np.empty(len(time_steps))

    for i, time in enumerate(time_steps):
        print(f"Timestep: {i}")
        mats = all_mats[i]

        blanket_mat = get_material_by_name(mats, material_name)

        dose_rates[i] = extract_contact_dose_rate(blanket_mat)

    output_variable = dose_rates

dose_rate_dict = {}

if __name__ == "__main__":
    if not use_stored_data:
        os.chdir(f"{base_dir}")

        U_dose_rates, Th_dose_rates = np.zeros((len(masses), num_timesteps)), np.zeros((len(masses), num_timesteps))

        for i, mass in enumerate(masses):
            print(f"~~~~~~ MASS: {mass} ~~~~~~~")

            os.chdir(f"Uranium/{mass}/{cooldow_folder_name}")
            U_results = Results("depletion_results.h5")
            time_steps = U_results.get_times()

            U_mats = []
            for j, time in enumerate(time_steps):
                U_mats.append(U_results.export_to_materials(j))
            os.chdir("../../..")

            os.chdir(f"Thorium/{mass}/{cooldow_folder_name}")
            Th_results = Results("depletion_results.h5")

            Th_mats = []
            for j, time in enumerate(time_steps):
                Th_mats.append(Th_results.export_to_materials(j))
            os.chdir("../../..")

            

            U_blanket_dose_rates, Th_blanket_dose_rates = 0, 0
            U_channel_dose_rates, Th_channel_dose_rates = 0, 0

            U_blanket_proc = Process(target=parallel_extract_dose_rate, 
                                    args=(U_mats, "doped flibe blanket", time_steps, U_blanket_dose_rates))
            U_channel_proc = Process(target=parallel_extract_dose_rate, 
                                   args=(U_mats, "doped flibe channels", time_steps, U_channel_dose_rates))
            
            Th_blanket_proc = Process(target=parallel_extract_dose_rate, 
                                    args=(Th_mats, "doped flibe blanket", time_steps, Th_blanket_dose_rates))
            Th_channel_proc = Process(target=parallel_extract_dose_rate, 
                                    args=(Th_mats, "doped flibe channels", time_steps, Th_channel_dose_rates))

            # Start all processes
            U_blanket_proc.start()
            U_channel_proc.start()
            Th_blanket_proc.start()
            Th_channel_proc.start()

            # Wait for all processes to finish before continuing
            U_blanket_proc.join()
            U_channel_proc.join()
            Th_blanket_proc.join()
            Th_channel_proc.join()

            U_dose_rates[i] = U_blanket_dose_rates + U_channel_dose_rates
            Th_dose_rates[i] = Th_blanket_dose_rates + Th_channel_dose_rates

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
        with open("dose_rate_dict.pkl", 'wb') as file:
            dose_rate_dict = pickle.load(file)


        try:
            os.chdir("figures")
        except:
            os.mkdir("figures")
            os.chdir("figures")

        for dopant in dopants:
            dose_rates = dose_rate_dict[dopant]
            time_steps = dose_rate_dict["time_steps"]

            fig, ax = plt.subplots()

            for i, dose_rate in enumerate(dose_rates):
                ax.plot(time_steps, dose_rates[i])

            ax.set_xlabel("Time (days)")
            ax.set_ylabel("Contact Dose Rate (Sv/hr)")
            ax.set_title(dopant)

            ax.set_xscale("log")

            fig.savefig(f"{dopant}_contact_dose_rate_decay.png", dpi=300)

            os.chdir("../..")



