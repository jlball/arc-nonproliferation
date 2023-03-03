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

def generate_device(dopant, dopant_mass, Li6_enrichment=7.5):
    device = anp.Device()

    # ==============================================================================
    # Geometry
    # ==============================================================================

    """ PFCs and Vacuum Vessel """

    vv_points = np.loadtxt("/home/jlball/arc-nonproliferation/data/arc_vv.txt")

    pfc_polygon = openmc.model.Polygon(vv_points, basis='rz')
    vv_inner_edge = pfc_polygon.offset(0.3) #PFC
    vv_channel_inner = vv_inner_edge.offset(1.0) #VV
    channel_outer = vv_channel_inner.offset(2.0) #FLiBe channels
    vv_channel_outer = channel_outer.offset(3.0) #Channel shell

    """ Blanket and Outer Blanket Tank """

    blanket_points = np.loadtxt("/home/jlball/arc-nonproliferation/data/arc_blanket.txt")

    blanket_inner = openmc.model.Polygon(blanket_points, basis='rz')
    blanket_outer = blanket_inner.offset(2) #Blanket tank outer

    regions = openmc.model.subdivide([pfc_polygon,
                                    vv_inner_edge, vv_channel_inner,
                                    channel_outer, vv_channel_outer,
                                    blanket_inner, blanket_outer])

    plasma, pfc, vv, channel, tank_inner, salt, tank_outer, outside = regions

    doped_flibe = anp.doped_flibe(dopant, dopant_mass, volume=1e8, Li6_enrichment=Li6_enrichment)

    device.plasma = openmc.Cell(region=plasma, fill=None, name='plasma')
    device.pfc = openmc.Cell(region=pfc, fill=anp.tungsten, name='PFC')
    device.vv = openmc.Cell(region=vv, fill=anp.vcrti_VV, name='VV')
    device.channel = openmc.Cell(region=channel, fill=doped_flibe, name='channels')
    device.tank_inner = openmc.Cell(region=tank_inner, fill=anp.vcrti_BI, name='tank inner')
    device.blanket = openmc.Cell(region=salt, fill=doped_flibe, name='blanket')
    device.tank_outer = openmc.Cell(region=tank_outer, fill=anp.vcrti_BO, name='tank outer')
    device.domain.region = device.domain.region & outside

    # ==============================================================================
    # Settings
    # ==============================================================================

    """ Source Definition """
    source = openmc.Source()
    source.space = openmc.stats.CylindricalIndependent(openmc.stats.Discrete(400, 1), openmc.stats.Uniform(a=-np.pi/18, b=np.pi/18), openmc.stats.Discrete(0, 1))
    source.angles = openmc.stats.Isotropic()
    source.energy = openmc.stats.Discrete([14.1E6], [1.0])

    device.settings.source = source
    device.settings.photon_transport = False

    device.settings.batches = 10

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
    flibe_filter = openmc.MaterialFilter(doped_flibe)
    device.add_tally('FLiBe Tally', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[flibe_filter])

    device.add_tally('Li Tally', ['(n,Xt)'], filters=[flibe_filter], nuclides=['Li6', 'Li7'])

    """ Breeding Tally """
    if dopant == 'U':
        device.add_tally('Breeding Tally', ['kappa-fission', 'absorption'], nuclides=['U238'], filters=[flibe_filter])
    if dopant == 'Th':
        device.add_tally('Breeding Tally', ['kappa-fission', 'absorption'], nuclides=['Th232'], filters=[flibe_filter])

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
device = generate_device('U', 0)
device.build()
plots.export_to_xml()
openmc.plot_geometry()

os.chdir("../..")

# ==============================================================================
# Scan
# ==============================================================================
particles = int(1e3)

Li6_enrichments = np.linspace(0, 100, num=5)
np.savetxt(base_dir + '/enrichments.txt', Li6_enrichments)

mass = 15e3

for enrichment in Li6_enrichments:
    U_device = generate_device('U', mass, Li6_enrichment=enrichment)
    Th_device = generate_device('Th', mass, Li6_enrichment=enrichment)

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



