import openmc
import openmc.deplete
import numpy as np
import os

###
# MATERIALS
###

flibe = openmc.Material()
flibe.add_elements_from_formula('F4Li2Be')
flibe.set_density('g/cm3', 1.94)
flibe.depletable = True

uf4 = openmc.Material()
uf4.add_elements_from_formula('UF4')
uf4.set_density('g/cm3', 6.7)
uf4.depletable = True

weight_percent = 0.01
doped_flibe = openmc.Material.mix_materials([uf4, flibe], [weight_percent, 1 - weight_percent], 'wo', name="doped flibe")
doped_flibe.volume = 4/3 * np.pi * np.power(100, 3)
doped_flibe.depletable = True

materials = openmc.Materials([doped_flibe])

###
# GEOMETRY
###

sphere_surface = openmc.Sphere(r=100, boundary_type='transmission')
bounding_sphere = openmc.Sphere(r=120, boundary_type='vacuum')

flibe_region = -sphere_surface
vacuum_region = +sphere_surface & -bounding_sphere

flibe_cell = openmc.Cell(fill=doped_flibe, region = flibe_region)
vacuum_cell = openmc.Cell(fill=None, region=vacuum_region)

geometry = openmc.Geometry([flibe_cell, vacuum_cell])

###
# SETTINGS
###
source = openmc.Source()
source.space = openmc.stats.Point()
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete(14.1e6, 1)

settings = openmc.Settings()
settings.photon_transport = False
settings.particles = int(1e3)
settings.batches = 10
settings.source = source
settings.run_mode = 'fixed source'

independent_device = openmc.model.Model(materials=materials, geometry=geometry, settings=settings)
coupled_device = openmc.model.Model(materials=materials, geometry=geometry, settings=settings)

## COUPLED DEPLETION ###
#Make a directory for and run the coupled depletion calculation

###
# DEPLETION SETTINGS
###
num_steps = 10
time_steps = [100*24*60*60 / num_steps] * num_steps
source_rates = [1e20] * num_steps

chain_file = '/home/jlball/arc-nonproliferation/data/chain_endfb71_pwr.xml'

# try:
#     os.mkdir("COUPLED")
#     os.chdir("COUPLED")
# except:
#     os.chdir("COUPLED")

# os.chdir('..')

# coupled_device.deplete(time_steps, 
#     source_rates=source_rates, 
#     operator_kwargs={'chain_file':chain_file, 
#                         'normalization_mode':'source-rate', 
#                         'reduce_chain':False,
#                         'reduce_chain_level':None}, 
#     directory="COUPLED")

### INDEPENDENT DEPLETION ###
try:
    os.mkdir("INDEPENDENT")
    os.chdir("INDEPENDENT")
except:
    os.chdir("INDEPENDENT")

independent_device.export_to_xml()

flux, micro_xs = openmc.deplete.get_microxs_and_flux(independent_device,
                                                        openmc.Materials.from_xml(),
                                                        chain_file=chain_file,
                                                        run_kwargs = {"threads":20,
                                                                    "particles":int(1e4)})

operator = openmc.deplete.IndependentOperator(openmc.Materials.from_xml(), 
                                                    flux,
                                                    micro_xs,
                                                    chain_file=chain_file, 
                                                    normalization_mode='source-rate', 
                                                    reduce_chain=False, 
                                                    reduce_chain_level=None, 
                                                    )

integrator = openmc.deplete.PredictorIntegrator(operator, 
                                                time_steps, 
                                                source_rates=source_rates,
                                                timestep_units='s')

integrator.integrate()
