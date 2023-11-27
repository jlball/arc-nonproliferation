import pickle
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.stats import linregress
import matplotlib.colors as colors
import matplotlib.cm as cm
import openmc
import openmc.plotter as plotter


# ====================================================
# Data loading
# ====================================================
masses = np.array([5, 10, 20, 30, 40, 50])
Li6_enrichments = np.array([2.5, 5, 7.5, 15, 30, 60, 90])
Li6_enrichments_str = np.array(['2.5', '5', '7.5', '15', '30', '60', '90'])
num_steps = 20

folder_prefix = 'pub_run_'

dopants = ["U", "Th"]

figure_folder = 'pub_figures'

try:
    os.mkdir(figure_folder)
except:
    print("skipping folder creation")

for dopant in dopants:
    time_to_sq = np.empty((len(masses), len(Li6_enrichments)))
    fission_power_t_sq = np.empty((len(masses), len(Li6_enrichments)))
    isotopic_purity = np.empty((len(masses), len(Li6_enrichments)))
    tbr_t_0 = np.empty((len(masses), len(Li6_enrichments)))
    flux_spectrum = np.empty((len(Li6_enrichments), len(masses), num_steps, 709, 2))
    reaction_spectra = np.empty((len(Li6_enrichments), len(masses), num_steps, 709, 2, 2))

    for i, enrichment in enumerate(Li6_enrichments_str):
        with open(folder_prefix + enrichment + f'/data/{dopant}_data_dict.pkl', 'rb') as file:
            data_dict = pickle.load(file)

            time_to_sq[:, i] = data_dict["time_to_sq"]/24
            fission_power_t_sq[:, i] = data_dict["fission_power_t_sq"]
            isotopic_purity[:, i] = data_dict["isotopic_purities"]
            tbr_t_0[:, i] = data_dict["tbr_t0"]
            flux_spectrum[i] = data_dict["flux_spectrum"]
            reaction_spectra[i] = data_dict["reaction_spectra"]


# ====================================================
# Plotting
# ====================================================
    try:
        os.chdir(figure_folder + f"/{dopant}")
    except:
        os.mkdir(figure_folder + f"/{dopant}")
        os.chdir(figure_folder + f"/{dopant}")

    if dopant == 'U':
        plt_color ='r'
        plt_cm = cm.Reds
    else:
        plt_color = 'g'
        plt_cm = cm.Greens

    norm = colors.Normalize(vmin=-0.5*masses[-1], vmax=1.1*masses[-1])
    enrichment_norm = colors.Normalize(vmin=-0.5*Li6_enrichments[-1], vmax=1.1*Li6_enrichments[-1])

    width_in = 7
    height_in = 7

    # +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
    # Time to 1 SQ

    fig, ax = plt.subplots()
    fig.set_size_inches(width_in,height_in)
    ax.spines["top"].set_color("None")
    ax.spines["right"].set_color("None")

    fontdict = {"size":18}

    for i, t_SQ in enumerate(time_to_sq):
        ax.scatter(Li6_enrichments, t_SQ, color=plt_cm(norm(masses[i])), s=15)
        ax.plot(Li6_enrichments, t_SQ, color=plt_cm(norm(masses[i])), alpha=0.5)
        #fit = linregress(Li6_enrichments, t_SQ)
        #ax.plot(Li6_enrichments, fit.slope*Li6_enrichments + fit.intercept, color=plt_color, alpha=0.3)

        #text_offset = (10, -4)
        #ax.annotate(f"{masses[i]} Tons", (Li6_enrichments[-1], t_SQ[-1]), color=plt_cm(norm(masses[i])), textcoords='offset points', xytext=text_offset)

    ax.set_yscale("log")
    ax.set_ylim(7, 3000)

    ax.set_xlabel("Li-6 Enrichment (percent)", fontdict=fontdict)
    ax.set_ylabel("Time to 1 SQ (days)", fontdict=fontdict)
    ax.set_title("Time to 1 SQ", fontdict=fontdict)
    
    fig.savefig(f"{dopant}_time_to_sq.png", dpi=300)
    fig.savefig(f"{dopant}_time_to_sq.svg")

    # +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
    # TBR

    fig, ax = plt.subplots()
    fig.set_size_inches(width_in,height_in)
    ax.spines["top"].set_color("None")
    ax.spines["right"].set_color("None")

    for i, tbr in enumerate(tbr_t_0):
        ax.scatter(Li6_enrichments, tbr, color=plt_cm(norm(masses[i])), s=15)
        ax.plot(Li6_enrichments, tbr, color=plt_cm(norm(masses[i])), alpha=0.5)

        #text_offset = 5
        #ax.annotate(f"{masses[i]} Tons", (Li6_enrichments[-1], tbr[-1]), color=plt_cm(norm(masses[i])), textcoords='offset points', xytext=(text_offset, 0))

    #ax.set_yscale("log")

    ax.set_xlabel("Li-6 Enrichment (percent)", fontdict=fontdict)
    ax.set_ylabel("TBR", fontdict=fontdict)
    ax.set_title("Tritium Breeding Ratio", fontdict=fontdict)

    ax.set_ylim(0.9, 1.25)

    fig.savefig(f"{dopant}_TBR.png", dpi=300)
    fig.savefig(f"{dopant}_TBR.svg")

    # +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
    # Isotopic Purity

    fig, ax = plt.subplots()
    fig.set_size_inches(width_in,height_in)
    ax.spines["top"].set_color("None")
    ax.spines["right"].set_color("None")

    for i, purity in enumerate(isotopic_purity):
        ax.scatter(Li6_enrichments, purity*100, color=plt_cm(norm(masses[i])), s=20)
        ax.plot(Li6_enrichments, purity*100, color=plt_cm(norm(masses[i])), alpha=0.5)

        #fit = linregress(Li6_enrichments, t_SQ)
        #ax.plot(Li6_enrichments, fit.slope*Li6_enrichments + fit.intercept, color=plt_color, alpha=0.3)

        #text_offset = (10, -4)
        #ax.annotate(f"{masses[i]} Tons", (Li6_enrichments[-1], purity[-1]), color=plt_cm(norm(masses[i])), textcoords='offset points', xytext=text_offset)

    ax.set_xlabel("Li-6 Enrichment (percent)", fontdict=fontdict)
    ax.set_ylabel("Percent Fissile Isotope (percent)", fontdict=fontdict)
    ax.set_title("Isotopic Purity", fontdict=fontdict)

    ax.set_ylim(98.9, 100)

    fig.savefig(f"{dopant}_isotopic_purity.png", dpi=300)
    fig.savefig(f"{dopant}_isotopic_purity.svg")

    # +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
    # Fission Power

    fig, ax = plt.subplots()
    fig.set_size_inches(width_in,height_in)
    ax.spines["top"].set_color("None")
    ax.spines["right"].set_color("None")

    for i, power in enumerate(fission_power_t_sq):
        ax.scatter(Li6_enrichments, power, color=plt_cm(norm(masses[i])), s=20)
        ax.plot(Li6_enrichments, power, color=plt_cm(norm(masses[i])), alpha=0.5)

        #fit = linregress(Li6_enrichments, t_SQ)
        #ax.plot(Li6_enrichments, fit.slope*Li6_enrichments + fit.intercept, color=plt_color, alpha=0.3)

        #text_offset = (10, -4)
        #ax.annotate(f"{masses[i]} Tons", (Li6_enrichments[-1], power[-1]), color=plt_cm(norm(masses[i])), textcoords='offset points', xytext=text_offset)

    ax.set_xlabel("Li-6 Enrichment (percent)", fontdict=fontdict)
    ax.set_ylabel("Fission Power (MW)", fontdict=fontdict)
    ax.set_title("Fission Power at $t_{SQ}$", fontdict=fontdict)

    ax.set_ylim(0, 120)
    
    fig.savefig(f"{dopant}_fission_Power_t_sq.png", dpi=300)
    fig.savefig(f"{dopant}_fission_Power_t_sq.svg")

    # +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
    # Flux Spectrum

    fig, ax = plt.subplots()
    ax.spines["top"].set_color("None")
    ax.spines["right"].set_color("None")

    energy_groups = openmc.mgxs.EnergyGroups(openmc.mgxs.GROUP_STRUCTURES['CCFE-709'])
    energy_bins = energy_groups.group_edges
    flux_energies = 0.5*(energy_bins[1:] + energy_bins[:-1])

    ax_xs = ax.twinx()
    ax_xs.spines["top"].set_color("None")

    if dopant == "U":
        xs_energies, xs_data = plotter.calculate_cexs("U238", ["capture"])
    if dopant == "Th":
        xs_energies, xs_data = plotter.calculate_cexs("Th232", ["capture"])

    ax_xs.plot(xs_energies, xs_data[0], color="grey", alpha=0.4)
    ax_xs.set_yscale("log")
    ax_xs.set_ylabel("Cross Section (barns)")

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
    fig.savefig(f"{dopant}_flux_spectrum.svg")

    # +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
    # (n, gamma)

    fig, ax = plt.subplots()
    ax.spines["top"].set_color("None")
    ax.spines["right"].set_color("None")

    for i, enrichment in enumerate(Li6_enrichments):
        ax.step(flux_energies, reaction_spectra[i, 0, 0, :, 1, 1], label=f"{enrichment} %", color=plt_cm(enrichment_norm(enrichment)))

    ax.set_xlabel("Neutron Energy (eV)")
    ax.set_ylabel("Reaction Rate (arc. units)")
    ax.set_title("(n, gamma) reaction rate vs. neutron energy")

    ax.set_yscale("log")
    ax.set_xscale("log")

    ax.set_xlim(1e-2, 20e6)

    fig.savefig(f"{dopant}_ngamma_spectrum.png")

    # +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+
    # fission

    fig, ax = plt.subplots()
    ax.spines["top"].set_color("None")
    ax.spines["right"].set_color("None")

    for i, enrichment in enumerate(Li6_enrichments):
        ax.step(flux_energies, reaction_spectra[i, 0, 0, :, 1, 0], label=f"{enrichment} %", color=plt_cm(enrichment_norm(enrichment)))

    ax.set_xlabel("Neutron Energy (eV)")
    ax.set_ylabel("Reaction Rate (arc. units)")
    ax.set_title("Fission reaction rate vs. neutron energy")

    ax.set_yscale("log")
    ax.set_xscale("log")

    ax.set_xlim(1e-2, 20e6)

    fig.savefig(f"{dopant}_fission_spectrum.png")

    os.chdir("../..")