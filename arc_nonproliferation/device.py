import numpy as np
import openmc
from arc_nonproliferation.components import *
from arc_nonproliferation.materials import *
from arc_nonproliferation.constants import neutrons_per_MJ
import shutil

"""
Class for simulating an ARC-class reactor blanket doped with either
U-238 or Th-233 for the purpose of breeding fissile material.

Instances of the class used for simulations are generated using the 
generate_device function which allows for a device instance to be generated 
with a given Li-6 enrichment and fertile inventory.

This class also stores as a class variable its fusion power and total neutron 
source rate so that these values can be easily accessed by both depletion scripts
and post processing scripts

"""
class Device(openmc.model.Model):

    fusion_power = 500 # MW
    neutron_source_rate = fusion_power * neutrons_per_MJ

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self._cells = []
        self._components = []
        self._tallies = []

        self.materials = []

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

    def build(self):
        """
        Builds the model geometry from specifically named components.
        Must be called before model can be run.
        """
        comp_cells = [cell for comp in self._components for cell in comp.cells]
        comp_regs = [~comp for comp in self._components if comp.exclude]

        #Exclude areas which have been defined as being components from total sim domain
        new_domain_reg = self.domain.region & openmc.Intersection(comp_regs)
        self.domain.region = new_domain_reg

        all_cells = comp_cells + self._cells
        self.univ = openmc.Universe(cells=all_cells)

        self.geometry = openmc.Geometry(self.univ)
        self.tallies = openmc.Tallies(self._tallies)

        if self.materials is []:
            self.materials = openmc.Materials(self.geometry.get_all_materials().values())
        
        super().export_to_xml()

    def add_tally(self, name, scores, nuclides=None, filters=None):
        """
        Creates a tally from given kwargs and adds it to the tallies object

        Parameters
        ----------
        name : str
            Name to be used for the tally
        scores : str
            tally scores, e. g. flux, absorption, (n,Xt)
        nuclides : list of str
            nuclides to individually track
        filters: list of openmc.Filter
                openmc filters like openmc.CellFilter to include in the tally

        Returns
        -------
        openmc.Tally, the generated tally object
        """
        tally = openmc.Tally(name=name)
        
        if nuclides is not None:
            tally.nuclides = nuclides

        tally.filters = filters
        tally.scores = scores

        self._tallies.append(tally)

        return tally

    def get_cell(self, name):
        comp_cells = [cell for comp in self._components for cell in comp.cells]
        all_cells = comp_cells + self._cells

        for cell in all_cells:
            if cell.name == name:
                return cell
        
        print("WARNING: Cell with name:", name, "not found. returning None.")
        return None

def generate_device(dopant, dopant_mass, Li6_enrichment=7.5, vv_file='arc_vv.txt', blanket_file="arc_blanket.txt", dopant_mass_units="kg"):
    """
    Generates a device object with specified fertile inventory and Li-6 enrichment.

    Parameters
    ----------
    dopant : str
        "U" for U-238 -> Pu-239, "Th" for Th-232 -> U233
    dopant_mass : float
        mass of fertile material to dope into the blanket in units of kg
    Li6_enrichment : float
        The percent of the lithium which is Li-6 instead of Li-7.
    vv_file : str
        name of file in the API data folder which contains the set of points which define the vv boundary
    blanket_file : str
        name of file in the API data folder which contains the set of points which define the blanket tank boundary

    Returns
    -------
    arc_nonproliferation.Device, the generated device object.
    """
    
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
    flibe_volume = vol_calc_load.volumes[8].n
    channels_volume = vol_calc_load.volumes[5].n


    if dopant_mass_units == "kg":
        doped_flibe = make_doped_flibe(dopant, 
                                        dopant_mass, 
                                        volume=flibe_volume + channels_volume, 
                                        Li6_enrichment=Li6_enrichment, 
                                        name="doped flibe blanket")
        
    elif dopant_mass_units == 'wppm':
        doped_flibe = make_impure_flibe(dopant, 
                                        dopant_mass, 
                                        name="doped flibe blanket")

    else:
        raise ValueError("Invalid units specified for dopant mass")
        
    device.doped_flibe = doped_flibe

    doped_flibe_channels = doped_flibe.clone()
    doped_flibe_channels.volume = channels_volume
    doped_flibe_channels.name = "doped flibe channels"

    doped_flibe_blanket = doped_flibe.clone()
    doped_flibe_blanket.volume = flibe_volume
    doped_flibe_blanket.name = "doped flibe blanket"

    vcrti_VV = vcrti.clone()
    vcrti.name = "VV"

    vcrti_BI = vcrti.clone()
    vcrti.name = "tank inner"

    vcrti_BO = vcrti.clone()
    vcrti.name = "tank outer"

    device.plasma = openmc.Cell(region=plasma, fill=None, name='plasma')
    device.pfc = openmc.Cell(region=pfc, fill=tungsten, name='PFC')
    device.vv = openmc.Cell(region=vv, fill=vcrti_VV, name='VV')
    device.channel = openmc.Cell(region=channel, fill=doped_flibe_channels, name='channel')
    device.tank_inner = openmc.Cell(region=tank_inner, fill=vcrti_BI, name='tank inner')
    device.blanket = openmc.Cell(region=salt, fill=doped_flibe_blanket, name='blanket')
    device.tank_outer = openmc.Cell(region=tank_outer, fill=vcrti_BO, name='tank outer')
    device.domain.region = device.domain.region & outside

    # ==============================================================================
    # Settings
    # ==============================================================================

    """ Source Definition """
    source = openmc.Source()
    source.space = openmc.stats.CylindricalIndependent(openmc.stats.Discrete(400, 1), openmc.stats.Uniform(a=0, b=2*np.pi), openmc.stats.Discrete(0, 1))
    source.angles = openmc.stats.Isotropic()
    source.energy = openmc.stats.Discrete([14.1E6], [1.0])

    device.settings.source = source

    return device

def generate_simple_device(dopant, dopant_mass, Li6_enrichment=7.5, vv_file='arc_vv.txt', blanket_file="arc_blanket.txt"):
    device = Device()
    device.dopant = dopant

    # ==============================================================================
    # Geometry
    # ==============================================================================

    """ PFCs and Vacuum Vessel """

    vv_points = np.loadtxt("/home/jlball/arc-nonproliferation/data/" + vv_file)

    pfc_polygon = openmc.model.Polygon(vv_points, basis='rz')
    vv_inner_edge = pfc_polygon.offset(0.3) #PFC
    vv_outer_edge = vv_inner_edge.offset(3.0) #Channel shell

    """ Blanket and Outer Blanket Tank """

    blanket_points = np.loadtxt("/home/jlball/arc-nonproliferation/data/" + blanket_file)

    blanket_inner = openmc.model.Polygon(blanket_points, basis='rz')
    gap = blanket_inner.offset(1.0)
    blanket_outer = gap.offset(2.0) #Blanket tank outer

    regions = openmc.model.subdivide([pfc_polygon,
                                    vv_inner_edge, 
                                    vv_outer_edge,
                                    blanket_inner, 
                                    blanket_outer])

    plasma, pfc, vv, salt, tank_outer, outside = regions

    # Read volume calc file
    vol_calc_load = openmc.VolumeCalculation.from_hdf5('/home/jlball/arc-nonproliferation/data/arc-1_volumes.h5')
    flibe_volume = vol_calc_load.volumes[8].n

    doped_flibe = make_doped_flibe(dopant, 
                                           dopant_mass, 
                                           volume=flibe_volume, 
                                           Li6_enrichment=Li6_enrichment, 
                                           name="doped flibe blanket")
    
    device.doped_flibe = doped_flibe

    doped_flibe_blanket = doped_flibe.clone()
    doped_flibe_blanket.volume = flibe_volume
    doped_flibe_blanket.name = "doped flibe blanket"

    vcrti_VV = vcrti.clone()
    vcrti.name = "VV"

    vcrti_BO = vcrti.clone()
    vcrti.name = "tank outer"

    device.plasma = openmc.Cell(region=plasma, fill=None, name='plasma')
    device.pfc = openmc.Cell(region=pfc, fill=tungsten, name='PFC')
    device.vv = openmc.Cell(region=vv, fill=vcrti_VV, name='VV')
    device.blanket = openmc.Cell(region=salt, fill=doped_flibe_blanket, name='blanket')
    device.tank_outer = openmc.Cell(region=tank_outer, fill=vcrti_BO, name='tank outer')
    device.domain.region = device.domain.region & outside

    # ==============================================================================
    # Settings
    # ==============================================================================

    """ Source Definition """
    source = openmc.Source()
    source.space = openmc.stats.CylindricalIndependent(openmc.stats.Discrete(400, 1), openmc.stats.Uniform(a=0, b=2*np.pi), openmc.stats.Discrete(0, 1))
    source.angles = openmc.stats.Isotropic()
    source.energy = openmc.stats.Discrete([14.1E6], [1.0])

    device.settings.source = source

    return device