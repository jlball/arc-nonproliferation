import openmc
from scipy.constants import Avogadro
from arc_nonproliferation.constants import *

"""
This module contains all material definitions and useful functions
for generating FLiBe materials doped with a specific mass of fertile
material.

"""

""" TUNGSTEN """
tungsten = openmc.Material(name='W')
tungsten.add_element('O',5/1e6,percent_type='wo')
tungsten.add_element('N',5/1e6,percent_type='wo')
tungsten.add_element('C',5/1e6,percent_type='wo')
tungsten.add_element('Na',4/1e6,percent_type='wo')
tungsten.add_element('K',2.5/1e6,percent_type='wo')
tungsten.add_element('Al',3/1e6,percent_type='wo')
tungsten.add_element('Ca',0.5/1e6,percent_type='wo')
tungsten.add_element('Cr',0.5/1e6,percent_type='wo')
tungsten.add_element('Cu',0.5/1e6,percent_type='wo')
tungsten.add_element('Fe',5/1e6,percent_type='wo')

tungsten.add_element('W',1-(5+5+5+4+2.5+3+0.5+0.5+0.5+5)/1e6,percent_type='wo')
tungsten.set_density('g/cm3',19.3)
#tungsten.volume = 65636 #Wedge vol
tungsten.depletable = False

""" FLIBE """
flibe = openmc.Material(name='FLiBe')
flibe.add_elements_from_formula('F4Li2Be')
flibe.set_density('g/cm3', 1.94)

""" V-4Cr-4Ti """
vcrti = openmc.Material(name='V-4Cr-4Ti VV')
vcrti.depletable = False

vcrti.add_element('Cr',0.04,percent_type='wo')
vcrti.add_element('Ti',0.04,percent_type='wo')

vcrti.add_element('C',56/1e6,percent_type='wo')
vcrti.add_element('O',181/1e6,percent_type='wo')
vcrti.add_element('N',103/1e6,percent_type='wo')
vcrti.add_element('B',7/1e6,percent_type='wo')
vcrti.add_element('Na',17/1e6,percent_type='wo')
vcrti.add_element('Mg',0.5/1e6,percent_type='wo')
vcrti.add_element('Al',119/1e6,percent_type='wo')
vcrti.add_element('Si',280/1e6,percent_type='wo')
vcrti.add_element('Mn',0.5/1e6,percent_type='wo')
vcrti.add_element('Fe',80/1e6,percent_type='wo')
vcrti.add_element('Ni',13/1e6,percent_type='wo')
vcrti.add_element('Cu',4/1e6,percent_type='wo')
vcrti.add_element('V',1-0.04-0.04-(56+181+103+7+17+0.5+119+280+0.5+80+13+4)/1e6,percent_type='wo')
vcrti.set_density('g/cm3',6.05) #This density value is sus and needs a good source

""" Uranium tetrafluroide """
uf4 = openmc.Material()
uf4.add_elements_from_formula('UF4')
uf4.set_density('g/cm3', 6.7)

""" Thorium tetrafluoride """
thf4 = openmc.Material()
thf4.add_elements_from_formula('ThF4')
thf4.set_density('g/cm3', 6.3)

""" Beryllium """
beryllium = openmc.Material()
beryllium.add_element("Be", 1.0)
beryllium.set_density('g/cm3', 1.85)

""" Fluorine """
fluorine = openmc.Material()
fluorine.add_element("F", 1.0)
fluorine.set_density("g/cm3", 0.001696)

""" Lithium """
lithium = openmc.Material()
lithium.add_element("Li", 1.0)
lithium.set_density("g/cm3", 0.534)

""" Uranium """
uranium = openmc.Material()
uranium.add_element("U", 1.0)
uranium.set_density('g/cm3', 19.1)

def get_tetrafluoride_mass(mass, dopant):
    """
    Computes mass of tetrafluroide from a given mass of pure dopant

    Parameters
    ----------
    mass : float
        mass of fertile material in grams
    dopant : str
        "U" for U-238 -> Pu-239, "Th" for Th-232 -> U233

    Returns
    -------
    float, mass of actinide tetrafluoride containing 'mass' grams of fertile material
    """

    if dopant == 'U':
        moles = mass / openmc.data.atomic_mass('U238')
        tetrafluoride_mass = moles * UF4_molar_mass

    elif dopant == 'Th':
        moles = mass / openmc.data.atomic_mass('Th232')
        tetrafluoride_mass = moles * ThF4_molar_mass
    else:
        raise ValueError("Not a valid dopant type")

    return tetrafluoride_mass

def make_doped_flibe(dopant, dopant_mass, Li6_enrichment=7.4, name='doped_flibe', volume=None, dopant_mass_units="kg"):
    """
    Return openmc material doped with specified fertile material

    Parameters
    ----------
    dopant : str
        "U" for U-238 -> Pu-239, "Th" for Th-232 -> U233
    dopant_mass : float
        mass of fertile material in kilograms
    Li6_enrichment : float
        The percent of the lithium which is Li-6 instead of Li-7.
    name : str
        the name of the material returned
    volume : the volume of the material returned in cubic centimeters
    
    Returns
    -------
    openmc.Material, FLiBe material doped with specified amount of fertile material
    """
    dopant_mass = dopant_mass * 1000 #Dopant mass is provided in units of kg, so here convert to grams

    flibe = openmc.Material()
    flibe.add_elements_from_formula('F4Li2Be', enrichment_target='Li6', enrichment_type='ao', enrichment=Li6_enrichment)
    flibe.set_density('g/cm3', 1.94)
    flibe.depletable = True

    if dopant == 'U':
        tetrafluoride = uf4
    elif dopant == 'Th':
        tetrafluoride = thf4
    else:
        raise ValueError("Invalid dopant passed into blanket liquid function")

    if dopant_mass_units == "kg":
        if volume == None:
            raise ValueError("Volume of blanket specified as None")
        else:
            flibe_mass = flibe.density * volume
            tetrafluoride_mass = get_tetrafluoride_mass(dopant_mass, dopant)
            tetrafluoride_weight_percent = tetrafluoride_mass / (flibe_mass + tetrafluoride_mass)
            
    elif dopant_mass_units == "wppm":
        tetrafluoride_weight_percent = dopant_mass/1e6
    
    else:
        raise ValueError("Invalid units given for dopant mass argument")

    doped_mat = openmc.Material.mix_materials([tetrafluoride, flibe], [tetrafluoride_weight_percent, 1 - tetrafluoride_weight_percent], 'wo', name=name)
    doped_mat.volume = volume
    doped_mat.depletable = True
    return doped_mat

def make_impure_flibe(wppm, name='doped_flibe'):
    
    weight_fraction = wppm/1e6

    impure_Be = openmc.Material.mix_materials([beryllium, uranium], [1 - weight_fraction, weight_fraction], 'wo')

    Be_fluoride = openmc.Material.mix_materials([impure_Be, fluorine], [1/3, 2/3], "ao")
    Li_fluoride = openmc.Material.add_elements_from_formula("LiF")

    impure_flibe = openmc.Material.mix_materials([Be_fluoride, Li_fluoride], [1/3, 2/3], "ao", name=name)

    impure_flibe.set_density("g/cm3", 1.94)

    return impure_flibe