import openmc


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

    