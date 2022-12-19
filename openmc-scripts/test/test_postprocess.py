import openmc
import arc_nonproliferation as anp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

sp = openmc.StatePoint('221219/statepoint.10.h5')
mesh_tally = sp.get_tally(name='Mesh Tally')


fusion_power = 500 #MW
total_source_rate = fusion_power * anp.neutrons_per_MJ # neutrons/s

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
print('Fission Rate per SP:', anp.get_uvalue(flibe_tally, 'fission'))
print('fission power assuming 200 MeV per reaction:', anp.get_uvalue(flibe_tally, 'fission') * total_source_rate * 3.204e-17) #200 MeV in MJ
print('kappa fission:', anp.get_uvalue(flibe_tally, 'kappa-fission') * total_source_rate * anp.MJ_per_eV, "MW")
print("fission q prompt:", anp.get_uvalue(flibe_tally, 'fission-q-prompt') * total_source_rate * anp.MJ_per_eV, "MW")
print('fission q recoverable:', anp.get_uvalue(flibe_tally, 'fission-q-recoverable') * total_source_rate * anp.MJ_per_eV, "MW")
print('heating:', anp.get_uvalue(flibe_tally, 'heating') * total_source_rate * anp.MJ_per_eV, "MW")
print('heating local:', anp.get_uvalue(flibe_tally, 'heating-local') * total_source_rate * anp.MJ_per_eV, "MW")

