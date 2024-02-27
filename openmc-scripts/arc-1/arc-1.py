import arc_nonproliferation as anp
import openmc
import numpy as np
import os
import sys

# ==============================================================================
# Geometry
# ==============================================================================

device = anp.generate_device("U", 20)

# Plotting
plot = openmc.Plot()
plot.filename = 'geometry_plot'
plot.basis = 'xz'
plot.origin = (350, 0, 0)
plot.width = (700, 800)
plot.pixels = (plot.width[0]*10, plot.width[1]*10)
plot.color_by = 'cell'

plots = openmc.Plots([plot])
plots.export_to_xml()


# ==============================================================================
# Settings
# ==============================================================================

""" Source Definition """
source = openmc.Source()
source.space = openmc.stats.CylindricalIndependent(openmc.stats.Discrete(450, 1), openmc.stats.Uniform(a=-np.pi/18, b=np.pi/18), openmc.stats.Discrete(0, 1))
source.angles = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([14.1E6], [1.0])

device.settings.source = source

# ==============================================================================
# Tallies
# ==============================================================================
# """ Cylindrical Mesh Tally """
# mesh = openmc.CylindricalMesh()
# mesh.r_grid = np.linspace(25, 200, num=25)
# mesh.z_grid = np.linspace(-200, 200, num=50)
# mesh.phi_grid = np.array([0, (2 * np.pi)/(18 * 2)])
# mesh_filter = openmc.MeshFilter(mesh)

# device.add_tally('Mesh Tally', ['flux', '(n,Xt)', 'heating-local', 'absorption'], filters=[mesh_filter])

# """ FLiBe Tally """
# flibe_filter = openmc.MaterialFilter(doped_mat)
# device.add_tally('FLiBe Tally', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[flibe_filter])

# ==============================================================================
# Run
# ==============================================================================

device.settings.photon_transport = True

device.build()
device.export_to_xml(remove_surfs=True)
openmc.plot_geometry()

#device.run(particles=int(1e3))

try:
    if sys.argv[1] is not None:
        os.mkdir(str(sys.argv[1]))
        device.move_files(str(sys.argv[1]))
        print("OpenMC files moved to new directory:", str(sys.argv[1]))

except:
    print("No directory specified, using this one")
