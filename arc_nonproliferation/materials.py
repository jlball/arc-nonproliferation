import openmc
from scipy.constants import Avogadro
from arc_nonproliferation.constants import *

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
tungsten.volume = 65636 #Wedge vol
tungsten.depletable = True

""" FLIBE """
flibe = openmc.Material(name='FLiBe')
flibe.add_elements_from_formula('F4Li2Be')
flibe.set_density('g/cm3', 1.94)

""" V-4Cr-4Ti """
vcrti_VV = openmc.Material(name='V-4Cr-4Ti VV')
vcrti_VV.volume = 221783 #Wedge vol
vcrti_VV.depletable = True
vcrti_BI = openmc.Material(name='V-4Cr-4Ti Blanket inner')
vcrti_BO = openmc.Material(name='V-4Cr-4Ti Blanket outer')
for m in [vcrti_VV,vcrti_BI,vcrti_BO]:
    #m.add_element('V',0.92,percent_type='wo')
    m.add_element('Cr',0.04,percent_type='wo')
    m.add_element('Ti',0.04,percent_type='wo')

    m.add_element('C',56/1e6,percent_type='wo')
    m.add_element('O',181/1e6,percent_type='wo')
    m.add_element('N',103/1e6,percent_type='wo')
    m.add_element('B',7/1e6,percent_type='wo')
    m.add_element('Na',17/1e6,percent_type='wo')
    m.add_element('Mg',0.5/1e6,percent_type='wo')
    m.add_element('Al',119/1e6,percent_type='wo')
    m.add_element('Si',280/1e6,percent_type='wo')
    m.add_element('Mn',0.5/1e6,percent_type='wo')
    m.add_element('Fe',80/1e6,percent_type='wo')
    m.add_element('Ni',13/1e6,percent_type='wo')
    m.add_element('Cu',4/1e6,percent_type='wo')
    m.add_element('V',1-0.04-0.04-(56+181+103+7+17+0.5+119+280+0.5+80+13+4)/1e6,percent_type='wo')
    m.set_density('g/cm3',6.05) #This density value is sus and needs a good source

"""Uranium tetrafluroide"""
uf4 = openmc.Material()
uf4.add_elements_from_formula('UF4')
uf4.set_density('g/cm3', 6.7)

"""Thorium tetrafluoride"""
thf4 = openmc.Material()
thf4.add_elements_from_formula('ThF4')
thf4.set_density('g/cm3', 6.3)

def get_tetrafluoride_mass(mass, dopant):
    """Computes mass of tetrafluroide from a given mass of pure dopant"""

    if dopant == 'U':
        moles = mass / openmc.data.atomic_mass('U238') * Avogadro
        tetrafluoride_mass = moles * UF4_molar_mass

    elif dopant == 'Th':
        moles = mass / openmc.data.atomic_mass('Th232') * Avogadro
        tetrafluoride_mass = moles * ThF4_molar_mass
    else:
        raise ValueError("Not a valid dopant type")

    return tetrafluoride_mass

def doped_flibe(dopant, dopant_mass, Li6_enrichment=7.8, name='doped_flibe', volume=None):
    """Return openmc material doped with specified fertile material"""

    dopant_mass = dopant_mass * 1000 #Dopant mass is provided in units of kg

    flibe = openmc.Material()
    flibe.add_elements_from_formula('F4Li2Be', enrichment_target='Li6', enrichment_type='ao', enrichment=Li6_enrichment)
    flibe.set_density('g/cm3', 1.94)

    if dopant == 'U':
        tetrafluoride = uf4
    elif dopant == 'Th':
        tetrafluoride = thf4
    else:
        raise ValueError("Non-implemented dopant passed into blanket liquid function")

    if volume == None:
        raise ValueError("Volume of blanket specified as None")
    else:
        flibe_mass = flibe.density * volume
        tetrafluoride_mass = get_tetrafluoride_mass(dopant_mass, dopant)
        tetrafluoride_weight_percent = tetrafluoride_mass / flibe_mass * 100

    doped_mat = openmc.Material.mix_materials([tetrafluoride, flibe], [tetrafluoride_weight_percent, 100 - tetrafluoride_weight_percent], 'wo', name="doped flibe")
    return doped_mat


    