import arc_nonproliferation.device as anp
import arc_nonproliferation.constants as anp_constants
import openmc
import openmc.deplete
import numpy as np
import os
import sys



def setup_device(device):
    """ Run settings """
    device.settings.photon_transport = False
    device.settings.particles = int(5e3)
    device.settings.batches = 100
    #device.survival_biasing = True

    device.settings.run_mode = "eigenvalue"

    """ Cell Filter """
    blanket_cell = device.get_cell(name='blanket')
    channel_cell = device.get_cell(name='channel')
    cell_filter = openmc.CellFilter([channel_cell, blanket_cell])

    """ Energy Filter """
    energy_filter = openmc.EnergyFilter.from_group_structure("CCFE-709")

    """ Mesh Filter """
    # mesh = openmc.CylindricalMesh()
    # mesh.r_grid = np.linspace(200, 700, num=500)
    # mesh.z_grid = np.linspace(-400, 400, num=800)
    # mesh.phi_grid = np.array([0, 2 * np.pi])
 
    # mesh_filter = openmc.MeshFilter(mesh)

    """ FLiBe Tally """
    device.add_tally('FLiBe Tally', ['(n,Xt)', 'fission', 'kappa-fission', 'heating', 'heating-local'], filters=[cell_filter])
    device.add_tally('Channel Flux Tally', ['flux'], filters=[energy_filter, cell_filter])
    device.add_tally('Li Tally', ['(n,Xt)'], filters=[energy_filter, cell_filter], nuclides=['Li6', 'Li7'])

    """ Fertile Tally """
    if device.dopant == "U":
        fertile_nuclide = "U238"
    elif device.dopant == "Th":
        fertile_nuclide = "Th232"
    else:
        raise ValueError("Invalid Dopant Type!")

    device.add_tally("Fertile Tally", ['fission','(n,gamma)'], filters=[energy_filter, cell_filter], nuclides=[fertile_nuclide])

    """ Mesh Tally """
    #device.add_tally("Mesh Filter", ['fission','(n,gamma)'], filters=[mesh_filter], nuclides=[fertile_nuclide])

    return device


U_device = setup_device(anp.generate_device("U", 50e3, Li6_enrichment=7.5))

U_device.run()
