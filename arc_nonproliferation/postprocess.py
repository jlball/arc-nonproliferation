import openmc
import numpy as np
import arc_nonproliferation as anp
import matplotlib.pyplot as plt
import uncertainties

def get_RZ_cyl_mesh_data(tally, score, value='mean', volume_norm=True):
    """Parse out a set of 3 2D np arrays for easy plotting of 
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
    data = data.reshape((r_points-1, z_points-1))

    if volume_norm:
        volumes = mesh.volumes.reshape((r_points-1, z_points-1)).T # Need to transpose matrix of volumes to get correct R/Z orientation
        data = data/volumes

    r_mesh, z_mesh = np.meshgrid(mesh.r_grid, mesh.z_grid)

    return r_mesh, z_mesh, data

def plot_RZ_quantity(tally, score, title='Title', volume_norm=True, cmap='plasma', value='mean'):
    """Plots an RZ quantity from a 2D cylindrical mesh

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
    ax.set_aspect(1)
    pcolormesh = ax.pcolormesh(r_grid, z_grid, tbr.T, cmap=cmap)
    fig.colorbar(pcolormesh, ax=ax, label='TBR')

    return fig, ax

def get_uvalue(tally, score, value='mean', filters=[]):
    """Gets a value and its uncertainty from a tally (without a mesh filter)

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