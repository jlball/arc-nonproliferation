import pickle
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.stats import linregress

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
    # Time to 1 SQ

    fig, ax = plt.subplots()
    ax.spines["top"].set_color("None")
    ax.spines["right"].set_color("None")

    for i, t_SQ in enumerate(time_to_sq):
        ax.scatter(Li6_enrichments, t_SQ, color=plt_color)
        ax.plot(Li6_enrichments, t_SQ, color=plt_color, alpha=0.3)

        #fit = linregress(Li6_enrichments, t_SQ)
        #ax.plot(Li6_enrichments, fit.slope*Li6_enrichments + fit.intercept, color=plt_color, alpha=0.3)

        text_offset = (10, -4)
        ax.annotate(f"{masses[i]} Tons", (Li6_enrichments[-1], t_SQ[-1]), color=plt_color, textcoords='offset points', xytext=text_offset)

    ax.set_yscale("log")

    ax.set_xlabel("Li-6 Enrichment (percent)")
    ax.set_ylabel("Time to 1 SQ (days)")
    

    fig.savefig(f"{dopant}_time_to_sq.png", dpi=300)

    # +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
    # TBR

    fig, ax = plt.subplots()
    ax.spines["top"].set_color("None")
    ax.spines["right"].set_color("None")

    for i, tbr in enumerate(tbr_t_0):
        ax.scatter(Li6_enrichments, tbr, color=plt_color)
        ax.plot(Li6_enrichments, tbr, color=plt_color, alpha=0.3)

        text_offset = 5
        ax.annotate(f"{masses[i]} Tons", (masses[-1], tbr[-1]), color=plt_color, textcoords='offset points', xytext=(text_offset, 0))

    ax.set_xlabel("Li-6 Enrichment (percent)")
    ax.set_ylabel("TBR")

    fig.savefig(f"{dopant}_TBR.png", dpi=300)

    # +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
    # Isotopic Purity

    fig, ax = plt.subplots()
    ax.spines["top"].set_color("None")
    ax.spines["right"].set_color("None")

    for i, purity in enumerate(isotopic_purity):
        ax.scatter(Li6_enrichments, purity, color=plt_color, s=10)
        ax.plot(Li6_enrichments, purity, color=plt_color, alpha=0.3)

        #fit = linregress(Li6_enrichments, t_SQ)
        #ax.plot(Li6_enrichments, fit.slope*Li6_enrichments + fit.intercept, color=plt_color, alpha=0.3)

        text_offset = (10, -4)
        ax.annotate(f"{masses[i]} Tons", (Li6_enrichments[-1], purity[-1]), color=plt_color, textcoords='offset points', xytext=text_offset)

    ax.set_xlabel("Li-6 Enrichment (percent)")
    ax.set_ylabel("Percent Fissile Isotope (percent)")

    fig.savefig(f"{dopant}_isotopic_purity.png", dpi=300)


    # +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
    # Fission Power

    fig, ax = plt.subplots()
    ax.spines["top"].set_color("None")
    ax.spines["right"].set_color("None")

    for i, power in enumerate(fission_power_t_sq):
        ax.scatter(Li6_enrichments, power, color=plt_color, s=10)
        ax.plot(Li6_enrichments, power, color=plt_color, alpha=0.3)

        #fit = linregress(Li6_enrichments, t_SQ)
        #ax.plot(Li6_enrichments, fit.slope*Li6_enrichments + fit.intercept, color=plt_color, alpha=0.3)

        text_offset = (10, -4)
        ax.annotate(f"{masses[i]} Tons", (Li6_enrichments[-1], power[-1]), color=plt_color, textcoords='offset points', xytext=text_offset)

    ax.set_xlabel("Li-6 Enrichment (percent)")
    ax.set_ylabel("Fission Power (MW)")

    fig.savefig(f"{dopant}_fission_Power_t_sq.png", dpi=300)