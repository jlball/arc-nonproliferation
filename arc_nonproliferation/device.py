import numpy as np
import openmc
from arc_nonproliferation.components import *

class Device(openmc.model.Model):
    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.cells = []
        self.components = []

        self.tallies = []

        # Define simulation boundary
        self.boundary = openmc.Sphere(r=1000, boundary_type='vacuum')
        self.domain = openmc.Cell(region=-self.boundary, name='domain')

        self._init_settings()

    def __setattr__(self, name, value):
        # Override __setattr__ to add components to the internal list
        if isinstance(value, Component):
            self.components.append(value)
        elif isinstance(value, openmc.Cell):
            self.cells.append(value)
        super().__setattr__(name, value)

    def _init_settings(self):
        """Initializes a settings object and sets it to fixed source mode"""
        settings = openmc.Settings()
        settings.run_mode = 'fixed source'
        self.settings = settings

    def build(self, model_angles=[-10, 10]):
        """Builds the model geometry from specifically named components"""
        comp_cells = [cell for comp in self.components for cell in comp.cells]
        comp_regs = [~comp for comp in self.components if comp.exclude]

        #Exclude areas which have been defined as being components from total sim domain
        new_domain_reg = self.domain.region & openmc.Intersection(comp_regs)
        self.domain.region = new_domain_reg

        all_cells = comp_cells + self.cells
        self.univ = openmc.Universe(cells=all_cells)

        #Creates wedge model with specified azimuthal width
        neg_plane = openmc.YPlane(boundary_type="reflective").rotate((0, 0, model_angles[0]))
        pos_plane = openmc.YPlane(boundary_type="reflective").rotate((0, 0, model_angles[1]))
        wedge_cell = openmc.Cell(region=+neg_plane & -pos_plane & -self.boundary,name='wedge cell')
        wedge_cell.fill = self.univ
        self.geometry = openmc.Geometry([wedge_cell])

    def add_tally(self, name, scores, nuclides=None, filters=None):
        """Creates a tally from given kwargs and adds it to the tallies object"""
        tally = openmc.Tally(name=name)
        
        if nuclides is not None:
            tally.nuclides = nuclides

        tally.filters = filters
        tally.scores = scores

        self.tallies.append(tally)

        return tally

    def run(self, threads=20, batches=10, particles=1000):
        """Runs the model with specified parameters"""
        self.settings.particles = particles
        self.settings.batches = batches

        self.statepointfile = super().run(threads=threads)