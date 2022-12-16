import arc_nonproliferation as anp
import openmc
import numpy as np

device = anp.Device()


"""PFCs and Vacuum Vessel:"""

vv_points = np.array([
    [100, 100],
    [100, -100],
    [-100, -100],
    [-100, 100]
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
    [-200, -200],
    [-200, 200]
])

blanket_inner = openmc.model.Polygon.from_file('blanket.txt', basis='rz')
blanket_outer = blanket_inner.offset(2) #Blanket tank outer

regions = openmc.model.subdivide([pfc_polygon,
                                  vv_inner_edge, vv_channel_inner,
                                  channel_outer, vv_channel_outer,
                                  blanket_inner, blanket_outer])

pfc, vv, channel, tank_inner, salt, tank_outer, outside = regions

device.pfc = openmc.Cell(region=pfc, fill=anp.tungsten, name='PFC')
device.vv = openmc.Cell(region=vv, fill=anp.vcrti_VV, name='VV')
device.channel = openmc.Cell(region=channel, fill=anp.flibe, name='channels')
device.tank_inner = openmc.Cell(region=tank_inner, fill=anp.vcrti_BI, name='tank inner')
device.blanket = openmc.Cell(region=salt, fill=anp.flibe, name='blanket')
device.tank_outer = openmc.Cell(region=tank_outer, fill=anp.vcrti_BO, name='tank outer')
device.domain.region = device.domain.region & outside

