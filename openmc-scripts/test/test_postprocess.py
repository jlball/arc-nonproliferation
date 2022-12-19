import openmc
import arc_nonproliferation as anp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

sp = openmc.StatePoint('statepoint.10.h5')
mesh_tally = sp.get_tally(name='Mesh Tally')


# =================================================================================
# Mesh TBR
# =================================================================================
fig, ax = anp.plot_RZ_quantity(mesh_tally, '(n,Xt)', volume_norm=False, title='Mesh TBR')
fig.set_figwidth(5)
fig.set_figheight(8)
fig.savefig('figures/mesh_tbr.png')

# =================================================================================
# FLiBe
# =================================================================================
flibe_tally = sp.get_tally(name='FLiBe Tally')

print('======== FLIBE TALLY VALUES =======')
print('Total TBR:', anp.get_uvalue(flibe_tally, '(n,Xt)'))