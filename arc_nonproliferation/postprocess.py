import openmc
from openmc.deplete import Results
import numpy as np
import arc_nonproliferation as anp
import matplotlib.pyplot as plt
import uncertainties
from numpy.polynomial.polynomial import Polynomial
from scipy.optimize import curve_fit, root
from scipy.interpolate import interp1d
import os

"""
This module houses functions which are useful for analysis of OpenMC
data across multiple different types of simulations (transport, 
coupled depletion, independent depletion)

"""

def get_RZ_cyl_mesh_data(tally, score, value='mean', volume_norm=True):
    """
    Parse out a set of 3 2D np arrays for easy plotting of 
    2D R/Z cylindrical Mesh tallies.

    Parameters
    ----------
    tally : openmc.Tally
        the tally object loaded from a statepoint file with a cylindrical mesh
    score : str
        The score you want plotted
    value : str, optional
        The value you want extracted, defaults to mean
    volume_norm : boolean, optional
        Whether or not to divide the value of each mesh element by the element volume


    Returns
    -------
    3 2D numpy arrays, containing the R coordinates, Z coordinates, and tally values at each location on the mesh respectively
    """
    # Get cylindrical mesh from the tally
    mesh = tally.find_filter(openmc.MeshFilter).mesh

    r_points = len(mesh.r_grid)
    z_points = len(mesh.z_grid)

    # Reshape the data into a 2 axis numpy array of R/Z flux values
    slice = tally.get_slice(scores=[score])
    data = slice.get_reshaped_data(value=value)
    data = data.reshape((z_points-1, r_points-1))

    if volume_norm:
        volumes = mesh.volumes.reshape((r_points-1, z_points-1)).T # Need to transpose matrix of volumes to get correct R/Z orientation
        data = data/volumes

    r_mesh, z_mesh = np.meshgrid(mesh.r_grid, mesh.z_grid)

    return r_mesh, z_mesh, data

def plot_RZ_quantity(tally, score, title='Title', volume_norm=True, cmap='plasma', value='mean'):
    """
    Plots an RZ quantity from a 2D cylindrical mesh

    Parameters
    ----------
    tally : openmc.Tally
        the tally object loaded from a statepoint file with a cylindrical mesh
    score : str
        The score you want plotted
    value : str, optional
        The value you want extracted, defaults to mean
    volume_norm : boolean, optional
        Whether or not to divide the value of each mesh element by the element volume
    cmap : str, optional
        Matplotlib color map to use in the plot


    Returns
    -------
    matplotlib fig and ax objects
    """
    r_grid, z_grid, tbr = anp.get_RZ_cyl_mesh_data(tally, score, volume_norm=volume_norm, value=value)

    fig, ax = plt.subplots()
    ax.set_xlabel('R (cm)')
    ax.set_ylabel('Z (cm)')
    ax.set_title(title)
    pcolormesh = ax.pcolormesh(r_grid, z_grid, tbr, cmap=cmap)
    fig.colorbar(pcolormesh, ax=ax, label='TBR')

    return fig, ax

def get_uvalue(tally, score, value='mean', filters=[]):
    """
    Gets a value and its uncertainty from a tally (without a mesh filter)

    Parameters
    ----------
    tally : openmc.Tally
        the tally object loaded from a statepoint file with a cylindrical mesh
    score : str
        The score you want plotted
    value : str, optional
        The value you want extracted, defaults to mean
    filters : list of OpenMC Filter objects
        A list of filters from which to pull the values

    Returns
    -------
    ufloat object, the desired value and its uncertainty
    """
    value = tally.get_values(scores=[score], value=value, filters=filters)
    std_dev = tally.get_values(scores=[score], value='std_dev', filters=filters)

    u_val = uncertainties.ufloat(value, std_dev)
    return u_val

def get_material_by_name(materials, name):
    """
    Gets a material object from a materials list from a name
    
    Parameters
    ----------
    materials : openmc.Materials
        The list of materials to search
    name : str
        The name of the material to return

    Returns
    -------
    the material with corresponding name
    """
    for mat in materials:
        if mat.name == name:
            return mat

def get_masses_from_mats(nuclide, results, density=False):
    """
    Returns mass of either Pu-239 or U-233 at each timestep from a depletion 
    results object. 

    Parameters
    ----------
    nuclide : str
        nuclide whose mass is to be extracted
    results : openmc.Results
        depletion results object to analyze

    Returns
    -------
    numpy.array, each entry being fissile mass in kg at each timestep
    """

    time_steps = results.get_times()
    masses = np.empty(len(time_steps))

    """ Extract list of fissile masses for each depletion time step """
    for i in range(0, len(time_steps)):
        materials = results.export_to_materials(i)

        doped_flibe_blanket = get_material_by_name(materials, 'doped flibe blanket')
        doped_flibe_channels = get_material_by_name(materials, 'doped flibe channels') 

        if density:
            mass = doped_flibe_blanket.get_mass_density(nuclide=nuclide) + doped_flibe_channels.get_mass_density(nuclide=nuclide)
        else:
            mass = doped_flibe_blanket.get_mass(nuclide=nuclide) + doped_flibe_channels.get_mass(nuclide=nuclide)

        masses[i] = mass / 1000 # Convert from grams to kg

        del materials
    
    return masses


def extract_time_to_sq(dopant, results):
    """
    Computes the time at which 1 significant quantity of fissile material is
    present in the blanket.

    Parameters
    ----------
    dopant : str
        "U" for U-238 -> Pu-239, "Th" for Th-232 -> U233
    results : openmc.Results
        the depletion results file to analyse

    Returns
    -------
    int, the index of the timestep just after t_SQ. used for later interpolation of relevant quantities at t_SQ

    float, the time in hours at which 1 SQ of fissile material is present in the blanket
    """

    time_steps = results.get_times(time_units='h')
    
    if dopant == "U":
        fissile_masses = get_masses_from_mats("Pu239", results)

    if dopant == "Th":
        fissile_masses = get_masses_from_mats("U233", results)

    # Get timestep with fissile mass nearest 1 SQ
    idx = np.abs(fissile_masses - anp.sig_quantity).argmin() 

    if fissile_masses[idx] > anp.sig_quantity:
        time_to_sq = np.interp(anp.sig_quantity, [fissile_masses[idx-1], fissile_masses[idx]], [time_steps[idx-1], time_steps[idx]])
        
    else:
        time_to_sq = np.interp(anp.sig_quantity, [fissile_masses[idx], fissile_masses[idx+1]], [time_steps[idx], time_steps[idx+1]])

    return idx, time_to_sq
    

def extract_decay_heat(results):
    """
    Gets a value and its uncertainty from a tally (without a mesh filter)

    Parameters
    ----------
    results : openmc.Results
        the depletion results file to analyse

    Returns
    -------
    numpy.array, decay heat power in MW at each time step.
    """

    time_steps = results.get_times()
    decay_heats = np.empty(len(time_steps))

    for i in range(0, len(time_steps)):
        materials = results.export_to_materials(i)
        doped_flibe = get_material_by_name(materials, 'doped flibe') 
        decay_heats[i] = doped_flibe.get_decay_heat()

    return decay_heats / 1e6 #Convert from watts to MW

def extract_isotopic_purity(dopant, results):
    """
    Extracts the purity of fissile material bred in the blanket

    Parameters
    ----------
    dopant : str
        "U" for U-238 -> Pu-239, "Th" for Th-232 -> U233
    results : openmc.Results
        the depletion results file to analyse
    idx : int
        the index of the timestep just after t_SQ

    Returns
    -------
    numpy.array, the ratio of bred fissile isotope (Pu-239 or U-233) to all other isotopes
    of Pu or U at each time step.
    """

    materials = results.export_to_materials(-1)

    doped_flibe_channels = get_material_by_name(materials, 'doped flibe channels')
    doped_flibe_blanket = get_material_by_name(materials, 'doped flibe blanket') 

    """ ~~~~~~~ Uranium -> Plutonium ~~~~~~ """
    if dopant == 'U':
        Pu_nuclides = doped_flibe_blanket.get_nuclides(element="Pu")
        atoms = {}
        total_atoms = 0
        for nuclide in Pu_nuclides:
            try:
                times, num_atoms_channels = results.get_atoms(mat=doped_flibe_channels, nuc=nuclide)
                times, num_atoms_blanket = results.get_atoms(mat=doped_flibe_blanket, nuc=nuclide)
            except:
                num_atoms_blanket = 0
                num_atoms_channels = 0
            atoms[nuclide] = num_atoms_channels + num_atoms_blanket
            total_atoms = total_atoms + atoms[nuclide]
        return atoms['Pu239']/total_atoms
    
    """ ~~~~~~~ Thorium -> Uranium ~~~~~~ """
    # Returns both U233 content AND U232 content
    if dopant == 'Th':
        U_nuclides = doped_flibe_blanket.get_nuclides(element="U")
        atoms = {}
        total_atoms = 0
        for nuclide in U_nuclides:
            try:
                times, num_atoms_channels = results.get_atoms(mat=doped_flibe_channels, nuc=nuclide)
                times, num_atoms_blanket = results.get_atoms(mat=doped_flibe_blanket, nuc=nuclide)
            except:
                num_atoms_blanket = 0
                num_atoms_channels = 0
            atoms[nuclide] = num_atoms_channels + num_atoms_blanket
            total_atoms = total_atoms + atoms[nuclide]
        return atoms['U233']/total_atoms, atoms["U232"]/total_atoms
    
def extract_activity(results, nuclide):
    timesteps = results.get_times()

    activities = np.empty(len(timesteps))

    for i, step in enumerate(timesteps):
        materials = results.export_to_materials(i)

        doped_flibe_channels = get_material_by_name(materials, 'doped flibe channels')
        doped_flibe_blanket = get_material_by_name(materials, 'doped flibe blanket')

        try:
            channels_act = doped_flibe_channels.get_activity(by_nuclide=True)
            blanket_act = doped_flibe_blanket.get_activity(by_nuclide=True)

            total_act = channels_act[nuclide] + blanket_act[nuclide]
        except:
            total_act = 0

        activities[i] = total_act
    
    return activities


def get_element_mass(material, element):

    mass = 0
    if len(element) == 1:
        for nuc in material.nuclides:
            if nuc.name[0] == element and nuc.name[1:].isnumeric(): #Checks to make sure that we only test 1 letter elements
                mass += material.get_mass(nuclide=nuc.name)

    elif len(element) == 2 and not nuc.name[1:].isnumeric():
        for nuc in material.nuclides:
            if nuc.name[0:2] == element:
                mass += material.get_mass(nuclide=nuc.name)

    return mass

def extract_contact_dose_rate(material):
    # Data in this file retrieved from: https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/air.html 
    # It has units of cm^2/g
    air_mu_en = np.loadtxt("/home/jlball/arc-nonproliferation/data/air_muen.txt")
    air_mu_en_energies = air_mu_en[:, 0]
    air_mu_en = air_mu_en[:, 2]
    air_mu_en_interp = interp1d(air_mu_en_energies, air_mu_en)

    decay_photon_dist = material.get_decay_photon_energy(units="Bq")
    mu_material = np.zeros(decay_photon_dist.x.shape)

    elements = material.get_elements()
    for element in elements:
        photon_data = openmc.data.IncidentPhoton.from_hdf5(f"/home/jlball/xs_data/endfb80_hdf5/photon/{element}.h5")

        reactions = photon_data.reactions
        rxn_keys = photon_data.reactions.keys()
        total_xs_data = np.zeros(decay_photon_dist.x.shape)

        for rxn_key in rxn_keys:
            total_xs_data += reactions[rxn_key].xs(decay_photon_dist.x) * 1e-24 # Convert from barns to cm^2

        try:
            mu = total_xs_data / (openmc.data.atomic_weight(element) * anp.atomic_mass_unit_grams)
            element_mass = get_element_mass(material, element)
            mu_material += (element_mass / material.get_mass()) * mu
            #print(element_mass)
        except:
            continue

    C = 3.6e9 * (1.602e-19)
    dose = 0
    for i, energy in enumerate(decay_photon_dist.x):
        dose +=  C * (air_mu_en_interp(energy/1e6)/mu_material[i]) * ((decay_photon_dist.p[i] * (energy/1e6))/(material.get_mass()/1e3))

    return dose

# def mass_attenuation_coeff():

# def extract_surface_dose_rate(material):
