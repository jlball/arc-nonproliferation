import arc_nonproliferation as anp
import openmc
import numpy as np
import os
import sys

if sys.argv[1] is not None:
    base_dir = str(sys.argv[1])
    os.mkdir(base_dir)
else:
    raise ValueError("No base directory specified!")

os.mkdir(base_dir + '/Uranium')
os.mkdir(base_dir + '/Thorium')

def setup_device(device):
    device.settings.photon_transport = False
    device.settings.batches = 100
    device.settings.particles = int(5e4)
    device.settings.survival_biasing = True

    """ Cell Filter """
    blanket_cell = device.get_cell(name='blanket')
    channel_cell = device.get_cell(name='channel')
    cell_filter = openmc.CellFilter([channel_cell, blanket_cell])

    """ Energy Filter """
    energy_filter = openmc.EnergyFilter.from_group_structure("CCFE-709")

    """ Mesh Filter """
    mesh = openmc.CylindricalMesh()
    mesh.r_grid = np.linspace(200, 700, num=500)
    mesh.z_grid = np.linspace(-400, 400, num=800)
    mesh.phi_grid = np.array([0, 2 * np.pi])

    mesh_filter = openmc.MeshFilter(mesh)

    #device.add_tally('Mesh Tally', ['flux', '(n,Xt)', 'heating-local', 'absorption'], filters=[mesh_filter])

    #device.add_tally('Li Tally', ['(n,Xt)'], filters=[], nuclides=['Li6', 'Li7'])

    """ Fertile Tally """
    if device.dopant == "U":
        fertile_nuclide = "U238"
    elif device.dopant == "Th":
        fertile_nuclide = "Th232"
    else:
        raise ValueError("Invalid Dopant Type!")

    device.add_tally("Fertile Tally", ['fission','(n,gamma)'], filters=[energy_filter, cell_filter], nuclides=[fertile_nuclide])

    """ Mesh Tally """
    device.add_tally("Mesh Filter", ['fission','(n,gamma)'], filters=[mesh_filter], nuclides=[fertile_nuclide])

    return device

# Plotting
# plot = openmc.Plot()
# plot.filename = 'geometry_plot'
# plot.basis = 'xz'
# plot.origin = (450, 0, 0)
# plot.width = (600, 600)
# plot.pixels = (2000, 2000)
# plot.color_by = 'cell'

# plots = openmc.Plots([plot])

# os.mkdir(base_dir + "/geometry_plots")
# os.chdir(base_dir + "/geometry_plots")

# generte arbitrary device
# device = anp.generate_device('U', 0)
# device.build()
# plots.export_to_xml()
# openmc.plot_geometry()

os.chdir("../..")

# ==============================================================================
# Scan
# ==============================================================================

masses = np.array([20e3])
np.savetxt(base_dir + '/masses.txt', masses)

for mass in masses:
    U_device = anp.generate_device('U', mass, Li6_enrichment=7.5)
    Th_device = anp.generate_device('Th', mass, Li6_enrichment=7.5)

    setup_device(U_device)
    setup_device(Th_device)

    print('=================== MASS:' + str(mass) + 'kg ===================')
    os.mkdir(base_dir + '/Uranium/' + str(mass))
    os.chdir(base_dir + '/Uranium/' + str(mass))
    U_device.build()
    U_device.run()
    os.chdir('../../..')

    os.mkdir(base_dir + '/Thorium/' + str(mass))
    os.chdir(base_dir + '/Thorium/' + str(mass))
    Th_device.build()
    Th_device.run()
    os.chdir('../../..')



