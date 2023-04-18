import arc_nonproliferation as anp
import openmc
import openmc.deplete
import numpy as np
import os
import sys

""" Handle command line arguments and setup directory structure """
if sys.argv[1] is not None:
    base_dir = str(sys.argv[1])
    os.mkdir(base_dir)
    print("OpenMC output saved to:", base_dir)

    os.mkdir(base_dir + '/Uranium')
    os.mkdir(base_dir + '/Thorium')

# This function handles the simulation specific setup of each device object
def setup_device(device):
    """ Run settings """
    device.settings.photon_transport = False
    device.settings.particles = int(1e3)
    device.settings.batches = 5

    """ Cylindrical Mesh Tally """
    mesh = openmc.CylindricalMesh()
    mesh.r_grid = np.linspace(100, 700, num=50)
    mesh.z_grid = np.linspace(-300, 300, num=50)
    mesh.phi_grid = np.array([0, (2 * np.pi)/(18 * 2)])
    mesh_filter = openmc.MeshFilter(mesh)

    # """ Spectrum Mesh """
    # spec_mesh = openmc.CylindricalMesh()
    # spec_mesh.r_grid = np.linspace(100, 700, num=50)
    # spec_mesh.z_grid = np.linspace(-15, 15, num=1)
    # spec_mesh.phi_grid = np.array([0, (2 * np.pi)/(18 * 2)])
    # spec_mesh_filter = openmc.MeshFilter(mesh)

    """ Cell Filter """
    blanket_cell = device.get_cell(name='blanket')
    blanket_filter = openmc.CellFilter(blanket_cell)

    """ Energy Filter """
    energy_filter = openmc.EnergyFilter.from_group_structure("CCFE-709")

    device.add_tally('Mesh Tally', ['flux', '(n,Xt)', 'heating-local', 'absorption'], filters=[mesh_filter])

    """ FLiBe Tally """
    #flibe_filter = openmc.MaterialFilter(anp.get_material_by_name(device.materials, "doped_flibe"))
    device.add_tally('FLiBe Tally', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[])
    device.add_tally('Flux Tally', ['flux'], filters=[energy_filter, blanket_filter])
    device.add_tally('Li Tally', ['(n,Xt)'], filters=[], nuclides=['Li6', 'Li7'])

    return device

# Plotting
poloidal_plot = openmc.Plot()
poloidal_plot.filename = 'poloidal_plane'
poloidal_plot.basis = 'xz'
poloidal_plot.origin = (400, 0, 0)
poloidal_plot.width = (700, 700)
poloidal_plot.pixels = (2000, 2000)
poloidal_plot.color_by = 'cell'

toroidal_plot = openmc.Plot()
toroidal_plot.filename = 'toroidal_plane'
toroidal_plot.basis = 'xy'
toroidal_plot.origin = (0, 0, 0)
toroidal_plot.width = (1300, 1300)
toroidal_plot.pixels = (2000, 2000)
toroidal_plot.color_by = 'cell'

plots = openmc.Plots([poloidal_plot, toroidal_plot])

os.mkdir(base_dir + "/geometry_plots")
os.chdir(base_dir + "/geometry_plots")

# generte arbitrary device
device = anp.generate_device('U', 0)
device.build()
plots.export_to_xml()
openmc.plot_geometry()

os.chdir("../..")

# ==============================================================================
# Depletion Run
# ==============================================================================

mass = np.array([100e3]) #kg of fertile material
np.savetxt(base_dir + '/mass.txt', mass)

""" DEPLETION SETTINGS """
print("~~~~~~~~~~~~~~~~~~ FERTILE MASS: " + str(mass) + " kg ~~~~~~~~~~~~~~~~~~")

fusion_power = 500 #MW
num_steps = 4
time_steps = [100*24*60*60 / num_steps] * num_steps
source_rates = [fusion_power * anp.neutrons_per_MJ] * num_steps

chain_file = '/home/jlball/arc-nonproliferation/data/simple_chain_endfb71_pwr.xml'

""" Generate blankets doped to specified mass """

U_device = setup_device(anp.generate_device("U", mass[0]))
Th_device = setup_device(anp.generate_device("Th", mass[0]))

""" Run depletion calculation """

# This is necessary to ensure the tally.xml file gets written into the directory where the depletion calculation will run
os.mkdir(base_dir + '/Uranium/'+ str(mass))
os.chdir(base_dir + '/Uranium/'+ str(mass))
U_device.build()

os.chdir('../../..')

U_device.deplete(time_steps, 
    source_rates=source_rates, 
    operator_kwargs={'chain_file':chain_file, 'normalization_mode':'source-rate'}, 
    directory=base_dir + '/Uranium/'+ str(mass))

os.mkdir(base_dir + '/Thorium/'+ str(mass))
os.chdir(base_dir + '/Thorium/'+ str(mass))
Th_device.build()
os.chdir('../../..')

Th_device.deplete(time_steps, 
    source_rates=source_rates, 
    operator_kwargs={'chain_file':chain_file, 'normalization_mode':'source-rate'}, 
    directory=base_dir + '/Thorium/' + str(mass))
