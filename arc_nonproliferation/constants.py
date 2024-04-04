from scipy.constants import Avogadro

neutrons_per_MWy = 1.118e25 # = 1 MW*year/17.6 MeV
neutrons_per_MJ = 3.546e17 # 1MJ/17.6 Mev

MJ_per_eV = 1.60218e-25 

UF4_molar_mass = 314.02 #g/mol
ThF4_molar_mass = 309.03 #g/mol

Pu239_mass_in_kg = 3.9695545e-25 #kg mass of 1 Pu239 nucleus in kg, from wolfram alpha
U233_mass_in_kg = 3.8697142e-25 #kg mass of 1 U233 nucleus in kg, from wolfram alpha

sig_quantity = 8 # kg, the same for both Pu239 and U233. Defined by the IAEA

chain_file = '/home/jlball/arc-nonproliferation/data/chain_endfb71_pwr.xml'
#chain_file = '/home/jlball/arc-nonproliferation/data/chain_tendl2019_jeff33.xml'

Np239_half_life = 2.356 #days
Pa233_half_life = 26.98 #days

atomic_mass_unit_grams = 1.66054e-24 #grams/amu

# Plotting configuration
u_marker = 'o'
th_marker = 's'

u_color = 'tab:orange'
th_color = 'tab:purple'

dpi = 300

title_y = 1.05

fontdict = {"size":16}

