import pickle
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import numpy as np

# ====================================================
# Data Storage
# ====================================================
masses = np.array([5, 10, 20, 30, 40, 50])
Li6_enrichments = np.array([2.5, 5, 7.5, 15, 30, 60, 90])
Li6_enrichments_str = np.array(['2.5', '5', '7.5', '15', '30', '60', '90'])
num_steps = 19

folder_prefix = 'pub_run_'

# ====================================================
# Uranium
# ====================================================

dopants = ["U", "Th"]

for dopant in dopants:
    time_to_sq = np.empty((len(masses), len(Li6_enrichments)))
    fission_power_t_sq = np.empty((len(masses), len(Li6_enrichments)))
    fission_power_t_0 = np.empty((len(masses), len(Li6_enrichments)))
    isotopic_purity = np.empty((len(masses), len(Li6_enrichments)))
    tbr_t_0 = np.empty((len(masses), len(Li6_enrichments)))
    flux_spectrum = np.empty(len(Li6_enrichments), len(masses), num_steps, 709, 2)

    for i, enrichment in enumerate(Li6_enrichments_str):
        with open(folder_prefix + enrichment + f'/data/{dopant}_data_dict.pkl', 'rb') as file:
            data_dict = pickle.load(file)

            time_to_sq[:, i] = data_dict["time_to_sq"]/24
            fission_power_t_sq[:, i] = data_dict["fission_power_t_sq"]
            fission_power_t_0[:, i] = data_dict["fission_power_t_0"]
            isotopic_purity[:, i] = data_dict["isotopic_purities"]
            tbr_t_0[:, i] = data_dict["tbr_t0"]
            flux_spectrum[i] = data_dict["flux_spectrum"]


    fig, ax = plt.subplots()
    fig.set_size_inches(8, 5)

    ax.set_xscale("log")

    X, Y = np.meshgrid(Li6_enrichments, masses)

    #cf = ax.contourf(X, Y, time_to_sq, levels=np.logspace(0, 4, num=1000), norm=colors.LogNorm(vmin=time_to_sq.min(), vmax=time_to_sq.max()), cmap=cm.plasma)
    #fig.colorbar(cf, label="Time to SQ (days)")

    # Time to 1 SQ
    cs_t_sq = ax.contour(X, Y, time_to_sq, levels=[14, 30, 60, 90, 180, 365], colors="red")
    ax.clabel(cs_t_sq, inline=True, fontsize=10)

    # Fission Power
    cs_fission_power = ax.contour(X, Y, fission_power_t_0, levels=[10, 25, 50, 75, 100], colors="green")
    ax.clabel(cs_fission_power, inline=True, fontsize=10)

    # TBR
    cs_tbr = ax.contour(X, Y, tbr_t_0, colors="blue")
    ax.clabel(cs_tbr, inline=True, fontsize=10)

    cs_purity = ax.contour(X, Y, isotopic_purity, colors='purple')
    ax.clabel(cs_purity, inline=True, fontsize=10)

    ax.set_xlabel("Li6 Enrichment (percent)")
    ax.set_ylabel("Fertile Mass (Metric Tons)")
    ax.set_title(f"{dopant} Proliferation Plot")

    fig.savefig(f"{dopant}_proliferation_plot.png", dpi=300)