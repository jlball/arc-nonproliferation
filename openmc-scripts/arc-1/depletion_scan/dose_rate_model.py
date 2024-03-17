import openmc

#######################################################
# Geometry
#######################################################

material_thickness = 100

surface_to_tally = 30
tally_thickness = 1

# Planes:
left_plane = openmc.XPlane(-10, boundary_type='reflective')
right_plane = openmc.XPlane(10, boundary_type='reflective')
top_plane = openmc.YPlane(10, boundary_type='reflective')
bot_plane = openmc.YPlane(-10, boundary_type='reflective')

back_plane = openmc.ZPlane(-material_thickness, boundary_type='vacuum')
surface_plane = openmc.ZPlane(0)

tally_surface = openmc.ZPlane(surface_to_tally)
tally_back = openmc.ZPlane(surface_to_tally + tally_thickness, boundary_type='vacuum')

# Regions:
walls = +left_plane & -right_plane & -top_plane & + bot_plane
material_reg = walls & +back_plane & -surface_plane
gap_reg = walls & +surface_plane & -tally_surface
tally_reg = walls & +tally_surface & -tally_back

# Cells:
material_cell = openmc.Cell(name="mat_cell", region=material_reg)
gap_cell = openmc.Cell(name="gap_cell", region=gap_reg)
tally_cell = openmc.Cell(name="tally_cell", region=tally_reg)

#######################################################
# Settings
#######################################################

settings = openmc.Settings()
settings.photon_transport = True
settings.batches = 100
settings.particles = 1e4
settings.run_mode = 'fixed source'

#######################################################
# Tallies
#######################################################

dose_tally = openmc.Tally(name="dose_tally")

cell_filter = openmc.CellFilter([tally_cell])

# These two lines taken from Jonathan Shimwell's fusion-neutronics repo
# https://github.com/fusion-energy/neutronics-workshop/blob/reshuffle/tasks/task_11_CSG_shut_down_dose_tallies/1_cell_based_shut_down_dose_rate_example.py
energies, pSv_cm2 = openmc.data.dose_coefficients(particle="photon", geometry="AP")
dose_filter = openmc.EnergyFunctionFilter(
    energies, pSv_cm2, interpolation="cubic"  # interpolation method recommended by ICRP
)

dose_tally.filters = [dose_tally, cell_filter]

def generate_dose_rate_model(material):
    material_cell.fill = material

    geometry = openmc.Geometry([material_cell, gap_cell, tally_cell])

    # Setup photon source
    source = openmc.IndependentSource()
    source.angle = openmc.stats.Isotropic()
    source.energy = material.get_decay_photon_energy()
    source.space = openmc.stats.Box([-10, -10, -material_thickness], [10, 10, 0])
    source.particle = 'photon'

    settings.source = source

    # Setup tallies
    tallies = openmc.Tallies(dose_tally)

    model = openmc.model.Model(geometry=geometry, settings=settings, tallies=tallies)

    return model 