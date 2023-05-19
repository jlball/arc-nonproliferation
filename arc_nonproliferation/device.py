import numpy as np
import openmc
from arc_nonproliferation.components import *
from arc_nonproliferation.materials import *
import shutil

class Device(openmc.model.Model):
    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self._cells = []
        self._components = []
        self._tallies = []

        self.dopant = None
        self.dopant_mass = 0

        self.name = 'Untitled'

        # Define simulation boundary
        self.boundary = openmc.Sphere(r=1000, boundary_type='vacuum')
        self.domain = openmc.Cell(region=-self.boundary, name='domain')

        self._init_settings()

    def __setattr__(self, name, value):
        # Override __setattr__ to add components to the internal list
        if isinstance(value, Component):
            self._components.append(value)
        elif isinstance(value, openmc.Cell):
            self._cells.append(value)
        super().__setattr__(name, value)

    def _init_settings(self):
        """Initializes a settings object and sets it to fixed source mode"""
        settings = openmc.Settings()
        settings.run_mode = 'fixed source'
        self.settings = settings

    def build(self, model_angles=[0, 360]):
        """Builds the model geometry from specifically named components"""
        comp_cells = [cell for comp in self._components for cell in comp.cells]
        comp_regs = [~comp for comp in self._components if comp.exclude]

        #Exclude areas which have been defined as being components from total sim domain
        new_domain_reg = self.domain.region & openmc.Intersection(comp_regs)
        self.domain.region = new_domain_reg

        all_cells = comp_cells + self._cells
        self.univ = openmc.Universe(cells=all_cells)

        #Creates wedge model with specified azimuthal width
        #neg_plane = openmc.YPlane(boundary_type="reflective").rotate((0, 0, model_angles[0]))
        #pos_plane = openmc.YPlane(boundary_type="reflective").rotate((0, 0, model_angles[1]))
        #wedge_cell = openmc.Cell(region=+neg_plane & -pos_plane & -self.boundary,name='wedge cell')
        #wedge_cell.fill = self.univ
        self.geometry = openmc.Geometry(self.univ)
        self.tallies = openmc.Tallies(self._tallies)
        super().export_to_xml()

    def add_tally(self, name, scores, nuclides=None, filters=None):
        """Creates a tally from given kwargs and adds it to the tallies object"""
        tally = openmc.Tally(name=name)
        
        if nuclides is not None:
            tally.nuclides = nuclides

        tally.filters = filters
        tally.scores = scores

        self._tallies.append(tally)

        return tally

    def dope_blanket(self):
        if self.dopant is None:
            doped_flibe = flibe
        else:
            doped_flibe = doped_flibe(self.dopant, self.dopant_mass)

    def run(self, threads=20, batches=10, particles=1000, cwd=None):
        """Runs the model with specified parameters"""
        self.settings.particles = particles
        self.settings.batches = batches

        self.statepointfile = super().run(threads=threads, cwd=cwd)
        return self.statepointfile

    def get_cell(self, name):
        comp_cells = [cell for comp in self._components for cell in comp.cells]
        all_cells = comp_cells + self._cells

        for cell in all_cells:
            if cell.name == name:
                return cell
        
        print("WARNING: Cell with name:", name, "not found. returning None.")
        return None

def generate_device(dopant, dopant_mass, Li6_enrichment=7.5, vv_file='arc_vv.txt', blanket_file="arc_blanket.txt"):
    device = Device()
    device.dopant = dopant

    # ==============================================================================
    # Geometry
    # ==============================================================================

    """ PFCs and Vacuum Vessel """

    vv_points = np.loadtxt("/home/jlball/arc-nonproliferation/data/" + vv_file)

    pfc_polygon = openmc.model.Polygon(vv_points, basis='rz')
    vv_inner_edge = pfc_polygon.offset(0.3) #PFC
    vv_channel_inner = vv_inner_edge.offset(1.0) #VV
    channel_outer = vv_channel_inner.offset(2.0) #FLiBe channels
    vv_channel_outer = channel_outer.offset(3.0) #Channel shell

    """ Blanket and Outer Blanket Tank """

    blanket_points = np.loadtxt("/home/jlball/arc-nonproliferation/data/" + blanket_file)

    blanket_inner = openmc.model.Polygon(blanket_points, basis='rz')
    gap = blanket_inner.offset(1.0)
    blanket_outer = gap.offset(2.0) #Blanket tank outer

    regions = openmc.model.subdivide([pfc_polygon,
                                    vv_inner_edge, vv_channel_inner,
                                    channel_outer, vv_channel_outer,
                                    blanket_inner, blanket_outer])

    plasma, pfc, vv, channel, tank_inner, salt, tank_outer, outside = regions

    # Read volume calc file
    vol_calc_load = openmc.VolumeCalculation.from_hdf5('/home/jlball/arc-nonproliferation/data/arc-1_volumes.h5')
    flibe_volume = vol_calc_load.volumes[7].n

    doped_flibe = make_doped_flibe(dopant, dopant_mass, volume=flibe_volume, Li6_enrichment=Li6_enrichment)
    device.doped_flibe = doped_flibe

    device.plasma = openmc.Cell(region=plasma, fill=None, name='plasma')
    device.pfc = openmc.Cell(region=pfc, fill=tungsten, name='PFC')
    device.vv = openmc.Cell(region=vv, fill=vcrti_VV, name='VV')
    device.channel = openmc.Cell(region=channel, fill=doped_flibe, name='channels')
    device.tank_inner = openmc.Cell(region=tank_inner, fill=vcrti_BI, name='tank inner')
    device.blanket = openmc.Cell(region=salt, fill=doped_flibe, name='blanket')
    device.tank_outer = openmc.Cell(region=tank_outer, fill=vcrti_BO, name='tank outer')
    device.domain.region = device.domain.region & outside

    # ==============================================================================
    # Settings
    # ==============================================================================

    """ Source Definition """
    source = openmc.Source()
    source.space = openmc.stats.CylindricalIndependent(openmc.stats.Discrete(475, 1), openmc.stats.Uniform(a=-np.pi/18, b=np.pi/18), openmc.stats.Discrete(0, 1))
    source.angles = openmc.stats.Isotropic()
    source.energy = openmc.stats.Discrete([14.1E6], [1.0])

    device.settings.source = source

    return device