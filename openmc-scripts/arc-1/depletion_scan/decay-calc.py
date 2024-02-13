import openmc
from openmc.deplete import Results, IndependentOperator
import sys
import os
import numpy as np
import pickle
from arc_nonproliferation.postprocess import *
from arc_nonproliferation.device import *


if sys.argv[1] is not None:
    base_dir = './' + sys.argv[1]
else:
    raise ValueError("No base directory specified!")

masses = np.loadtxt(base_dir + "/masses.txt")

with open(f'{base_dir}/data/U_data_dict.pkl', 'rb') as file:
    data_dict = pickle.load(file)

    U_time_to_sq = data_dict["time_to_sq"]/24

#
#  Uranium
#

os.chdir(f"{base_dir}/Uranium")

for i, mass in enumerate(masses):
    os.chdir(f"{mass}")

    results = Results("depletion_results.h5")

    blanket_mat = lin_interp_material(results, "doped flibe blanket", U_time_to_sq[i])
    channel_mat = lin_interp_material(results, "doped flibe channels", U_time_to_sq[i])

    # Perform depletion calc

    device = generate_device("U", 0, Li6_enrichment=7.5)

    operator = IndependentOperator(openmc.Materials([device.channel.fill, device.blanket.fill]),
                                    fluxes=[0, 0],
                                    micros=[openmc.deplete.MicroXS(), openmc.deplete.MicroXS()],
                                    normalization_mode="source-rate")

    os.chdir("..")

