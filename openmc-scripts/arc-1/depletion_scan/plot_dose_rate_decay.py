import openmc
from openmc.deplete import Results
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import sys
import os
import pickle
from arc_nonproliferation.postprocess import get_material_by_name
from arc_nonproliferation.constants import *
from scipy.optimize import root
from Li6_plots import norm

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

sp_fig, sp_ax = plt.subplots()

# ====================================================
# Data extraction
# ====================================================

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
    time_steps = data_dict["dose_rate_cooldown_times"]

    #When this func is zero, dose rate is 1 Sv/hr, NRC self protecting limit
    def interp_log_dose_decay(time, mass_idx):
        if time > 0:
            return np.interp(time, time_steps, np.log10(dose_rate_decay[int(mass_idx)]))
        else:
            return np.log10(dose_rate_decay[int(mass_idx)][0]) - 100 * time

    self_protecting_times = np.zeros(len(masses))
    for i, mass in enumerate(masses):
        if np.max(dose_rate_decay[i]) < 1:
            self_protecting_times[i] = 0
        else:
            res = root(interp_log_dose_decay, np.array([10]), args=(i))
            self_protecting_times[i] = res.x

    data_dict["self_protecting_time"] = self_protecting_times

    # ====================================================
    # Plotting
    # ====================================================
    plot_masses = masses/1e3

    os.chdir(f"{base_dir}/figures")

    if dopant == 'U':
        plt_cm = cm.Oranges
    else:
        plt_cm = cm.Purples

    # Plot dose rate decay
    fig, ax = plt.subplots()
    ax.spines["top"].set_color("None")
    ax.spines["right"].set_color("None")

    for i, mass in enumerate(plot_masses):
        ax.loglog(activities[0]/(3600*24), dose_rate_decay[i], label=f"{mass}", color=plt_cm(norm(mass)))

    ax.hlines([1], activities[0][0]/(3600*24), activities[0][-1]/(3600*24), color="red", linestyles="dashed")

    ax.set_xlabel("Time (days)", fontdict=fontdict)
    ax.set_ylabel("Dose rate at 1 m (Sv/hr)", fontdict=fontdict)

    fig.savefig(f"{dopant}_dose_rate_decay.png", dpi=dpi)
    fig.savefig(f"{dopant}_dose_rate_decay.pdf")

    # Plot self protection time

    sp_ax.spines["top"].set_color("None")
    sp_ax.spines["right"].set_color("None")

    if dopant == "U":
        sp_ax.scatter(plot_masses, self_protecting_times, label="U-238", marker=u_marker, color=u_color)
    else:
        sp_ax.scatter(plot_masses, self_protecting_times, label="Th-232", marker=th_marker, color=th_color)

    sp_ax.set_ylabel("Self-protection Time (days)", fontdict=fontdict)
    sp_ax.set_xlabel("Fertile Inventory (metric tons)", fontdict=fontdict)

    
    # ====================================================
    # Save data
    # ====================================================

    os.chdir("../..")

    with open(f'{base_dir}/data/{dopant}_data_dict.pkl', 'wb') as file:
        pickle.dump(data_dict, file)

# Save self protection plot
os.chdir(f"{base_dir}/figures")
sp_fig.savefig("self_protecting_time.png", dpi=300)
sp_fig.savefig("self_protecting_time.pdf")

    

