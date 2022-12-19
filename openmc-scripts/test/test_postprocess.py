import openmc
import arc_nonproliferation as anp
import numpy as np

sp = openmc.StatePoint('statepoint.10.h5')
mesh_tally = sp.get_tally(name='Mesh Tally')

