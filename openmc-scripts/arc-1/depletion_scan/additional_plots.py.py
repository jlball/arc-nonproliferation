import pickle
import matplotlib.pyplot as plt
import numpy as np
import os

# ====================================================
# Data loading
# ====================================================
masses = np.array([5, 10, 20, 30, 40, 50])
Li6_enrichments = np.array([2.5, 5, 7.5, 15, 30, 60, 90])
Li6_enrichments_str = np.array(['2.5', '5', '7.5', '15', '30', '60', '90'])

folder_prefix = 'pub_run_'

dopants = ["U", "Th"]

for dopant in dopants:
    time_to_sq = np.empty((len(masses), len(Li6_enrichments)))
    fission_power_t_sq = np.empty((len(masses), len(Li6_enrichments)))
    isotopic_purity = np.empty((len(masses), len(Li6_enrichments)))
    tbr_t_0 = np.empty((len(masses), len(Li6_enrichments)))

    for i, enrichment in enumerate(Li6_enrichments_str):
        with open(folder_prefix + enrichment + f'/data/{dopant}_data_dict.pkl', 'rb') as file:
            U_data_dict = pickle.load(file)

            time_to_sq[:, i] = U_data_dict["time_to_sq"]/24
            fission_power_t_sq[:, i] = U_data_dict["fission_power_t_sq"]
            isotopic_purity[:, i] = U_data_dict["isotopic_purities"]
            tbr_t_0[:, i] = U_data_dict["tbr_t0"]


# ====================================================
# Plotting
# ====================================================

    if dopant == 'U':
        plt_color ='r'
    else:
        plt_color = 'g'


    # +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
    # TBR

    fig, ax = plt.subplots()
    ax.spines["top"].set_color("None")
    ax.spines["right"].set_color("None")

    for i, tbr in enumerate(tbr_t_0.T):
        ax.scatter(masses, tbr, color=plt_color)
        ax.plot(masses, tbr, color=plt_color)

        text_offset = 5
        ax.annotate(f"{Li6_enrichments_str[i]} %", (masses[-1], tbr[-1]), color=plt_color, textcoords='offset points', xytext=(text_offset, 0))

    ax.set_xlabel("Fertile Mass (Metric Tons)")
    ax.set_ylabel("TBR")

    fig.savefig(f"{dopant}_TBR.png", dpi=300)
