import openmc
from openmc.deplete import Results
import numpy as np
import arc_nonproliferation as anp
import matplotlib.pyplot as plt
import uncertainties
from numpy.polynomial.polynomial import Polynomial
from scipy.optimize import curve_fit, root
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

def get_masses_from_mats(dopant, results):
    """
    Returns mass of either Pu-239 or U-233 at each timestep from a depletion 
    results object. 

    Parameters
    ----------
    dopant : str
        "U" for uranium doped blanket, returns Pu-239 mass, or "Th" for a
        thorium doped blanket, returns U-233 mass
    results : openmc.Results
        depletion results object to analyze

    Returns
    -------
    numpy.array, each entry being fissile mass in kg at each timestep
    """

    time_steps = results.get_times()
    fissile_masses = np.empty(len(time_steps))

    """ Extract list of fissile masses for each depletion time step """
    for i in range(0, len(time_steps)):
        materials = results.export_to_materials(i)

        doped_flibe_blanket = get_material_by_name(materials, 'doped flibe blanket')
        doped_flibe_channels = get_material_by_name(materials, 'doped flibe channels') 

        if dopant == "U":
            fissile_mass = doped_flibe_blanket.get_mass(nuclide='Pu239') + doped_flibe_channels.get_mass(nuclide='Pu239')
        elif dopant == "Th":
            fissile_mass = doped_flibe_blanket.get_mass(nuclide='U233') + doped_flibe_channels.get_mass(nuclide='U233')
        else:
            raise ValueError("Invalid dopant type passed into extract time to SQ function")

        fissile_masses[i] = fissile_mass / 1000 # Convert from grams to kg
    
    return fissile_masses

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
    float, the time in hours at which 1 SQ of fissile material is present in the blanket
    """

    time_steps = results.get_times(time_units='h')
    
    fissile_masses = get_masses_from_mats(dopant, results)

    # Get timestep with fissile mass nearest 1 SQ
    idx = np.abs(fissile_masses - 8).argmin() 

    """ Linear fit to fissile masses data to determine time to SQ """
    fit = Polynomial.fit(time_steps[idx-1:idx + 1], fissile_masses[idx-1:idx + 1], 1)
    time_to_sig_quantity = (fit - anp.sig_quantity).roots()[0]
    return time_to_sig_quantity

def extract_time_to_sq_curve_fit(dopant, results):
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
    float, the time in hours at which 1 SQ of fissile material is present in the blanket
    """
    def fit_func(x, A, B, C):
        return A*np.exp(-B*x) + (C * x) - A

    time_steps = results.get_times(time_units='h')

    fissile_masses = get_masses_from_mats(dopant, results)

    # Get timestep with fissile mass nearest 1 SQ
    idx = np.abs(fissile_masses - 8).argmin() 


    p0 = np.array([1, 1, 1])
    """ fit to fissile masses data to determine time to SQ """
    popt, pcov = curve_fit(fit_func, time_steps, fissile_masses, p0)

    """ get root """

    res = root(fit_func, 1000, args=(popt[0], popt[1], popt[2]-8))

    return res.x
    

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

    Returns
    -------
    numpy.array, the ratio of bred fissile isotope (Pu-239 or U-233) to all other isotopes
    of Pu or U at each time step.
    """

    materials = results.export_to_materials(-1)
    doped_flibe_channels = get_material_by_name(materials, 'doped flibe channels')
    doped_flibe_blanket = get_material_by_name(materials, 'doped flibe blanket') 

    """ ~~~~~~~ Uranium -> Plutonium ~~~~~~ """
    Pu_nuclides = ["Pu238", "Pu239", "Pu240", "Pu241"]
    if dopant == 'U':
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
    U_nuclides =["U232", "U233", "U234", "U235", "U236", "U237", "U238"]
    if dopant == 'Th':
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
        return atoms['U233']/total_atoms
    