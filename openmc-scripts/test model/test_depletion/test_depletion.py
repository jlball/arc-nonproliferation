import arc_nonproliferation as anp
import openmc
import openmc.deplete
import numpy as np
import sys
import os

if sys.argv[1] is not None:
    base_dir = './' + sys.argv[1]
else:
    raise ValueError("No base directory specified!")

os.mkdir(base_dir)
os.chdir(base_dir)

device = anp.Device()

# ==============================================================================
# Geometry
# ==============================================================================

""" PFCs and Vacuum Vessel """

vv_points = np.array([
    [100, 100],
    [100, -100],
    [50, -100],
    [50, 100]
])

pfc_polygon = openmc.model.Polygon(vv_points, basis='rz')
vv_inner_edge = pfc_polygon.offset(0.3) #PFC
vv_channel_inner = vv_inner_edge.offset(1.0) #VV
channel_outer = vv_channel_inner.offset(2.0) #FLiBe channels
vv_channel_outer = channel_outer.offset(3.0) #Channel shell

""" Blanket and Outer Blanket Tank """

blanket_points = np.array([
    [200, 200],
    [200, -200],
    [25, -200],
    [25, 200]
])

blanket_inner = openmc.model.Polygon(blanket_points, basis='rz')
blanket_outer = blanket_inner.offset(2) #Blanket tank outer

regions = openmc.model.subdivide([pfc_polygon,
                                  vv_inner_edge, vv_channel_inner,
                                  channel_outer, vv_channel_outer,
                                  blanket_inner, blanket_outer])

plasma, pfc, vv, channel, tank_inner, salt, tank_outer, outside = regions

doped_mat = anp.doped_flibe('U', 1e4, volume=1e8)

device.plasma = openmc.Cell(region=plasma, fill=None, name='plasma')
device.pfc = openmc.Cell(region=pfc, fill=anp.tungsten, name='PFC')
device.vv = openmc.Cell(region=vv, fill=anp.vcrti_VV, name='VV')
device.channel = openmc.Cell(region=channel, fill=doped_mat, name='channels')
device.tank_inner = openmc.Cell(region=tank_inner, fill=anp.vcrti_BI, name='tank inner')
device.blanket = openmc.Cell(region=salt, fill=doped_mat, name='blanket')
device.tank_outer = openmc.Cell(region=tank_outer, fill=anp.vcrti_BO, name='tank outer')
device.domain.region = device.domain.region & outside

# ==============================================================================
# Settings
# ==============================================================================

""" Source Definition """
source = openmc.Source()
source.space = openmc.stats.CylindricalIndependent(openmc.stats.Discrete(75, 1), openmc.stats.Uniform(a=-np.pi/18, b=np.pi/18), openmc.stats.Discrete(0, 1))
source.angles = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([14.1E6], [1.0])

device.settings.source = source

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
flibe_filter = openmc.MaterialFilter(doped_mat)
device.add_tally('FLiBe Tally', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[flibe_filter])

# ==============================================================================
# Depletion Run
# ==============================================================================

""" DEPLETION SETTINGS """
fusion_power = 500 #MW
src_rate = fusion_power * anp.neutrons_per_MJ

num_steps = 10
time_steps = [365*34*60*60 / num_steps] * num_steps
source_rates = [src_rate, src_rate, src_rate, 0, 0, 0, 0, 0, 0, 0] 

""" Setup device model """
device.settings.photon_transport = False
device.settings.particles = int(1e3)
device.settings.batches = 5

device.build()
device.export_to_xml(remove_surfs=True)

""" Setup depletion solver """
operator = openmc.deplete.CoupledOperator(device, '../chain_endfb71_pwr.xml', normalization_mode='source-rate')
integrator = openmc.deplete.PredictorIntegrator(operator=operator, timesteps=time_steps, source_rates=source_rates)

""" Run depletion calculation """
integrator.integrate()