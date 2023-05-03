import openmc
import arc_nonproliferation as anp
import os
import shutil

dir = "geometry_plots"

if os.path.exists(dir):
    shutil.rmtree(dir)

os.mkdir(dir)
os.chdir(dir)

# generte arbitrary device
device = anp.generate_device('U', 0)
device.build()

print(device.materials)

# Plotting
poloidal_plot = openmc.Plot()
poloidal_plot.filename = 'poloidal_plane'
poloidal_plot.basis = 'xz'
poloidal_plot.origin = (400, 0, 0)
poloidal_plot.width = (700, 700)
poloidal_plot.pixels = (2000, 2000)
poloidal_plot.color_by = 'material'
poloidal_plot.colors = {
    anp.tungsten: (36, 36, 36),
    anp.vcrti_VV: (140, 140, 140),
    anp.vcrti_BI: (140, 140, 140),
    anp.vcrti_BO: (140, 140, 140),
    device.doped_flibe: (82, 255, 209),
}

toroidal_plot = openmc.Plot()
toroidal_plot.filename = 'toroidal_plane'
toroidal_plot.basis = 'xy'
toroidal_plot.origin = (0, 0, 0)
toroidal_plot.width = (1300, 1300)
toroidal_plot.pixels = (2000, 2000)
toroidal_plot.color_by = 'cell'

plots = openmc.Plots([poloidal_plot, toroidal_plot])

plots.export_to_xml()
openmc.plot_geometry()