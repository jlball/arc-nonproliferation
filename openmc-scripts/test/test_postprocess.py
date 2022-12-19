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
r_grid, z_grid, tbr = anp.get_RZ_cyl_mesh_data(mesh_tally, '(n,Xt)', volume_norm=False)

fig, ax = plt.subplots()
ax.set_xlabel('R (cm)')
ax.set_ylabel('Z (cm)')
ax.set_title("TBR Mesh Plot")
pcolormesh = ax.pcolormesh(r_grid, z_grid, tbr.T)
fig.colorbar(pcolormesh, ax=ax, label='TBR')
ax.set_aspect(1)
plt.savefig('figures/mesh_tbr.png')