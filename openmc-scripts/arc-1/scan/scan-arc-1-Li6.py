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

    # ==============================================================================
    # Tallies
    # ==============================================================================
    """ Cylindrical Mesh Tally """
    mesh = openmc.CylindricalMesh()
    mesh.r_grid = np.linspace(25, 200, num=25)
    mesh.z_grid = np.linspace(-200, 200, num=50)
    mesh.phi_grid = np.array([0, (2 * np.pi)/(18 * 2)])
    mesh_filter = openmc.MeshFilter(mesh)

    device.add_tally('Mesh Tally', ['flux', '(n,Xt)', 'heating-local', 'absorption'], filters=[mesh_filter])

    """ FLiBe Tally """
    device.add_tally('FLiBe Tally', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[])

    device.add_tally('Li Tally', ['(n,Xt)'], filters=[], nuclides=['Li6', 'Li7'])

    """ Breeding Tally """
    print(device.dopant)
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
particles = int(1e3)

#Li6_enrichments = np.linspace(0, 100, num=5)
Li6_enrichments = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 50.0, 70.0, 100.0]
#Li6_enrichments = [7.5, 20.0]
np.savetxt(base_dir + '/enrichments.txt', Li6_enrichments)

mass = 15e3

for enrichment in Li6_enrichments:
    U_device = setup_device(anp.generate_device('U', mass, Li6_enrichment=enrichment))
    Th_device = setup_device(anp.generate_device('Th', mass, Li6_enrichment=enrichment))

    print('=================== ENRICHMENT:' + str(enrichment) + ' ===================')
    os.mkdir(base_dir + '/Uranium/' + str(enrichment))
    os.chdir(base_dir + '/Uranium/' + str(enrichment))
    U_device.build()
    U_device.run(particles=particles)
    os.chdir('../../..')

    os.mkdir(base_dir + '/Thorium/' + str(enrichment))
    os.chdir(base_dir + '/Thorium/' + str(enrichment))
    Th_device.build()
    Th_device.run(particles=particles)
    os.chdir('../../..')



