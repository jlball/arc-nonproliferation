import openmc
from openmc.deplete import Results
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import pickle
from arc_nonproliferation.postprocess import get_material_by_name

cooldow_folder_name = "dose_rate_cooldown"

#openmc.config['chain_file'] = chain_file

dopants = ["U", "Th"]

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

        activities = results.get_activity(get_material_by_name(mats, "doped flibe blanket"))

        ax.loglog(activities[0]/(3600*24), activities[1], label=f"{mass}")
        os.chdir("../..")

    fig.savefig("dose_rate_act_decay.png", dpi=300)
    os.chdir("../..")