import arc_nonproliferation as anp
import openmc
import numpy as np

device = anp.Device()

"""PFCs and Vacuum Vessel:"""

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

"""Blanket and Outer Blanket Tank:"""

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

"""Source Definition"""
source = openmc.Source()
source.space = openmc.stats.CylindricalIndependent(openmc.stats.Discrete(75, 1), openmc.stats.Uniform(a=-np.pi/18, b=np.pi/18), openmc.stats.Discrete(0, 1))
source.angles = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([14.1E6], [1.0])

device.settings.source = source

device.build()
device.export_to_xml(remove_surfs=True)

device.run()