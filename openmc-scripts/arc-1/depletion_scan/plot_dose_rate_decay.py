import openmc
from openmc.deplete import Results
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import pickle
from arc_nonproliferation.postprocess import get_material_by_name
from scipy.optimize import root

cooldow_folder_name = "dose_rate_cooldown"
folder_prefix = 'pub_data_Li6_'

#openmc.config['chain_file'] = chain_file

dopants = ["U", "Th"]
enrichment = 7.5

decay_time_steps = np.logspace(-1, 4, num=8)

if sys.argv[1] is not None:
    base_dir = './' + sys.argv[1]
else:
    raise ValueError("No base directory specified!")

masses = np.loadtxt(base_dir + "/masses.txt")

for dopant in dopants:
    if dopant == "U":
        os.chdir(f"{base_dir}/Uranium")
    elif dopant == "Th":
        os.chdir(f"{base_dir}/Thorium")
    else:
        raise ValueError("Invalid dopant type")

    fig, ax = plt.subplots()

    for i, mass in enumerate(masses):
        os.chdir(f"{mass}")
        os.chdir(f"{cooldow_folder_name}")

        results = Results("depletion_results.h5")
        mats = results.export_to_materials(0)

        activities = results.get_activity(get_material_by_name(mats, "doped flibe blanket"), units="Bq")

        ax.loglog(activities[0]/(3600*24), activities[1], label=f"{mass}")
        os.chdir("../..")

    os.chdir("../..")
    os.chdir(f"{base_dir}/figures")
    fig.savefig(f"{dopant}_dose_rate_act_decay.png", dpi=300)
    os.chdir("../..")

    with open(f'{base_dir}/data/{dopant}_data_dict.pkl', 'rb') as file:
        data_dict = pickle.load(file)

    dose_rate_decay = data_dict["dose_rate_cooldown"]
    #dose_rate_act = data_dict["dose_rate_cooldown_act"]
    #self_protecting_time = data_dict["self_protecting_time"]
    time_steps = data_dict["dose_rate_cooldown_times"]

    os.chdir(f"{base_dir}/figures")

    # Plot dose rate decay
    fig, ax = plt.subplots()

    for i, mass in enumerate(masses):
        ax.loglog(activities[0]/(3600*24), dose_rate_decay[i], label=f"{mass}")

    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Dose rate at 30cm (Sv/hr)")

    ax.legend()
    fig.savefig(f"{dopant}_dose_rate_decay.png", dpi=300)
    
    
    #When this func is zero, dose rate is 1 Sv/hr, NRC self protecting limit
    def interp_log_dose_decay(time, mass_idx):
        return np.interp(time, time_steps, np.log10(dose_rate_decay[int(mass_idx)]))

    self_protecting_times = np.zeros(len(masses))
    for i, mass in enumerate(masses):
        res = root(interp_log_dose_decay, np.array([300]), args=(i))
        self_protecting_times[i] = res.x

    data_dict["self_protecting_time"] = self_protecting_times

    fig, ax = plt.subplots()

    ax.scatter(masses/1e3, self_protecting_times)
    ax.set_ylabel("Self-protecting Time (days)")
    ax.set_xlabel("Fertile Inventory (metric tons)")

    fig.savefig(f"{dopant}_self_protecting_time.png", dpi=300)

    os.chdir("../..")

    # # Plot activity decay
    # fig, ax = plt.subplots()

    # for i, mass in enumerate(masses):
    #     ax.loglog(activities[0]/(3600*24), dose_rate_act[i], label=f"{mass}")

    # ax.set_title("Total activity")
    # ax.set_ylabel("Time (days)")

    # ax.legend()

    # os.chdir(f"{base_dir}/figures")
    # fig.savefig(f"{dopant}_dose_rate_act_decay_from_calc.png", dpi=300)
    # os.chdir("../..")

    with open(f'{base_dir}/data/{dopant}_data_dict.pkl', 'wb') as file:
        pickle.dump(data_dict, file)

    

