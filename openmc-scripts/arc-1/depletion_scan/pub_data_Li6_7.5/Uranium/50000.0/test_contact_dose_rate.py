import openmc
import openmc.deplete
import arc_nonproliferation as anp


openmc.config['chain_file'] = anp.chain_file

res = openmc.deplete.Results("depletion_results.h5")

mats = res.export_to_materials(-1)

print(anp.surface_dose_rate(mats[1]))