import arc_nonproliferation as anp
import openmc
import numpy as np
import os
import sys

if sys.argv[1] is not None:
    base_dir = './' + sys.argv[1]
else:
    raise ValueError("No base directory specified!")

os.mkdir(base_dir)
#os.chdir(base_dir)

os.mkdir(base_dir + '/Uranium')
os.mkdir(base_dir + '/Thorium')

def setup_device(device):
    device.settings.photon_transport = False
    device.settings.batches = 10
    device.settings.particles = int(1e5)

    """ Cylindrical Mesh Tally """
    mesh = openmc.CylindricalMesh()
    mesh.r_grid = np.linspace(25, 200, num=25)
    mesh.z_grid = np.linspace(-200, 200, num=50)
    mesh.phi_grid = np.array([0, (2 * np.pi)/(18 * 2)])
    mesh_filter = openmc.MeshFilter(mesh)

    device.add_tally('Mesh Tally', ['flux', '(n,Xt)', 'heating-local', 'absorption'], filters=[mesh_filter])

    device.add_tally('Li Tally', ['(n,Xt)'], filters=[], nuclides=['Li6', 'Li7'])

    """ Breeding Tally """
    if device.dopant == 'U':
        device.add_tally('Breeding Tally', ['kappa-fission', 'absorption'], nuclides=['U238'], filters=[])
    if device.dopant == 'Th':
        device.add_tally('Breeding Tally', ['kappa-fission', 'absorption'], nuclides=['Th232'], filters=[])

    return device

# Plotting
plot = openmc.Plot()
plot.filename = 'geometry_plot'
plot.basis = 'xz'
plot.origin = (450, 0, 0)
plot.width = (600, 600)
plot.pixels = (2000, 2000)
plot.color_by = 'cell'

plots = openmc.Plots([plot])

os.mkdir(base_dir + "/geometry_plots")
os.chdir(base_dir + "/geometry_plots")

# generte arbitrary device
device = anp.generate_device('U', 0)
device.build()
plots.export_to_xml()
openmc.plot_geometry()

os.chdir("../..")

# ==============================================================================
# Scan
# ==============================================================================

masses = np.array([5e3, 10e3])
np.savetxt(base_dir + '/masses.txt', masses)

for mass in masses:
    U_device = anp.generate_device('U', mass)
    Th_device = anp.generate_device('Th', mass)

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



