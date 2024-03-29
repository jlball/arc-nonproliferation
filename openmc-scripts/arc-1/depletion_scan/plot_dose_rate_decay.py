import openmc
from openmc.deplete import Results
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import pickle
from arc_nonproliferation.postprocess import get_material_by_name

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
    dose_rate_act = data_dict["dose_rate_cooldown_act"]

    # Plot dose rate decay
    fig, ax = plt.subplots()

    for i, mass in enumerate(masses):
        ax.loglog(activities[0]/(3600*24), dose_rate_decay[i], label=f"{mass}")

    ax.legend()

    os.chdir(f"{base_dir}/figures")
    fig.savefig(f"{dopant}_dose_rate_decay.png", dpi=300)
    os.chdir("../..")

    # Plot activity decay
    fig, ax = plt.subplots()

    for i, mass in enumerate(masses):
        ax.loglog(activities[0]/(3600*24), dose_rate_act[i], label=f"{mass}")

    ax.set_title("Total activity")
    ax.set_ylabel("Time (days)")

    ax.legend()

    os.chdir(f"{base_dir}/figures")
    fig.savefig(f"{dopant}_dose_rate_act_decay_from_calc.png", dpi=300)
    os.chdir("../..")

    

