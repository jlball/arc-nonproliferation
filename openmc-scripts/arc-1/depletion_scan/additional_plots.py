import pickle
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.stats import linregress
import matplotlib.colors as colors
import matplotlib.cm as cm
import openmc


# ====================================================
# Data loading
# ====================================================
masses = np.array([5, 10, 20, 30, 40, 50])
Li6_enrichments = np.array([2.5, 5, 7.5, 15, 30, 60, 90])
Li6_enrichments_str = np.array(['2.5', '5', '7.5', '15', '30', '60', '90'])
num_steps = 20

folder_prefix = 'pub_run_'

dopants = ["U", "Th"]

for dopant in dopants:
    time_to_sq = np.empty((len(masses), len(Li6_enrichments)))
    fission_power_t_sq = np.empty((len(masses), len(Li6_enrichments)))
    isotopic_purity = np.empty((len(masses), len(Li6_enrichments)))
    tbr_t_0 = np.empty((len(masses), len(Li6_enrichments)))
    flux_spectrum = np.empty((len(Li6_enrichments), len(masses), num_steps, 709, 2))

    for i, enrichment in enumerate(Li6_enrichments_str):
        with open(folder_prefix + enrichment + f'/data/{dopant}_data_dict.pkl', 'rb') as file:
            data_dict = pickle.load(file)

            time_to_sq[:, i] = data_dict["time_to_sq"]/24
            fission_power_t_sq[:, i] = data_dict["fission_power_t_sq"]
            isotopic_purity[:, i] = data_dict["isotopic_purities"]
            tbr_t_0[:, i] = data_dict["tbr_t0"]
            flux_spectrum[i] = data_dict["flux_spectrum"]

# ====================================================
# Plotting
# ====================================================

    if dopant == 'U':
        plt_color ='r'
        plt_cm = cm.Reds
    else:
        plt_color = 'g'
        plt_cm = cm.Greens

    norm = colors.Normalize(vmin=-0.5*masses[-1], vmax=1.1*masses[-1])
    enrichment_norm = colors.Normalize(vmin=-0.5*Li6_enrichments[-1], vmax=1.1*Li6_enrichments[-1])

    # +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
    # Time to 1 SQ

    fig, ax = plt.subplots()
    ax.spines["top"].set_color("None")
    ax.spines["right"].set_color("None")

    for i, t_SQ in enumerate(time_to_sq):
        ax.scatter(Li6_enrichments, t_SQ, color=plt_cm(norm(masses[i])), s=15)
        ax.plot(Li6_enrichments, t_SQ, color=plt_cm(norm(masses[i])), alpha=0.5)

        #fit = linregress(Li6_enrichments, t_SQ)
        #ax.plot(Li6_enrichments, fit.slope*Li6_enrichments + fit.intercept, color=plt_color, alpha=0.3)

        text_offset = (10, -4)
        ax.annotate(f"{masses[i]} Tons", (Li6_enrichments[-1], t_SQ[-1]), color=plt_cm(norm(masses[i])), textcoords='offset points', xytext=text_offset)

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
        ax.scatter(Li6_enrichments, tbr, color=plt_cm(norm(masses[i])))
        ax.plot(Li6_enrichments, tbr, color=plt_cm(norm(masses[i])))

        text_offset = 5
        ax.annotate(f"{masses[i]} Tons", (Li6_enrichments[-1], tbr[-1]), color=plt_cm(norm(masses[i])), textcoords='offset points', xytext=(text_offset, 0))

    #ax.set_yscale("log")

    ax.set_xlabel("Li-6 Enrichment (percent)")
    ax.set_ylabel("TBR")

    fig.savefig(f"{dopant}_TBR.png", dpi=300)

    # +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
    # Isotopic Purity

    fig, ax = plt.subplots()
    ax.spines["top"].set_color("None")
    ax.spines["right"].set_color("None")

    for i, purity in enumerate(isotopic_purity):
        ax.scatter(Li6_enrichments, purity, color=plt_cm(norm(masses[i])), s=10)
        ax.plot(Li6_enrichments, purity, color=plt_cm(norm(masses[i])))

        #fit = linregress(Li6_enrichments, t_SQ)
        #ax.plot(Li6_enrichments, fit.slope*Li6_enrichments + fit.intercept, color=plt_color, alpha=0.3)

        text_offset = (10, -4)
        ax.annotate(f"{masses[i]} Tons", (Li6_enrichments[-1], purity[-1]), color=plt_cm(norm(masses[i])), textcoords='offset points', xytext=text_offset)

    ax.set_xlabel("Li-6 Enrichment (percent)")
    ax.set_ylabel("Percent Fissile Isotope (percent)")

    fig.savefig(f"{dopant}_isotopic_purity.png", dpi=300)


    # +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
    # Fission Power

    fig, ax = plt.subplots()
    ax.spines["top"].set_color("None")
    ax.spines["right"].set_color("None")

    for i, power in enumerate(fission_power_t_sq):
        ax.scatter(Li6_enrichments, power, color=plt_cm(norm(masses[i])), s=10)
        ax.plot(Li6_enrichments, power, color=plt_cm(norm(masses[i])), alpha=0.3)

        #fit = linregress(Li6_enrichments, t_SQ)
        #ax.plot(Li6_enrichments, fit.slope*Li6_enrichments + fit.intercept, color=plt_color, alpha=0.3)

        text_offset = (10, -4)
        ax.annotate(f"{masses[i]} Tons", (Li6_enrichments[-1], power[-1]), color=plt_cm(norm(masses[i])), textcoords='offset points', xytext=text_offset)

    ax.set_xlabel("Li-6 Enrichment (percent)")
    ax.set_ylabel("Fission Power (MW)")

    fig.savefig(f"{dopant}_fission_Power_t_sq.png", dpi=300)

    # +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
    # Flux Spectrum

    fig, ax = plt.subplots()
    ax.spines["top"].set_color("None")
    ax.spines["right"].set_color("None")

    energy_groups = openmc.mgxs.EnergyGroups(openmc.mgxs.GROUP_STRUCTURES['CCFE-709'])
    energy_bins = energy_groups.group_edges
    flux_energies = 0.5*(energy_bins[1:] + energy_bins[:-1])

    for i, enrichment in enumerate(Li6_enrichments):
        ax.step(flux_energies, flux_spectrum[i, 0, 0, :, 1], label=f"{enrichment} %", color=plt_cm(enrichment_norm(enrichment)))

    ax.set_xlabel("Neutron Energy (eV)")
    ax.set_ylabel("Average Flux (arb. units)")
    ax.set_title("Average Flux in BLanket Tank")

    ax.set_yscale("log")
    ax.set_xscale("log")

    ax.set_xlim(1e-2, 20e6)
    ax.set_ylim(1e-5, 1e1)

    fig.savefig(f"{dopant}_flux_spectrum.png", dpi=300)