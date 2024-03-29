import arc_nonproliferation as anp
import openmc
import openmc.deplete
import numpy as np
import os
import sys

""" Handle command line arguments and setup directory structure """
if sys.argv[1] is not None:
    base_dir = str(sys.argv[1])
    os.mkdir(base_dir)
    print("OpenMC output saved to:", base_dir)

    os.mkdir(base_dir + '/Uranium')
    os.mkdir(base_dir + '/Thorium')

if sys.argv[2] is not None:
    mass = float(sys.argv[2]) * 1e3

else:
    raise RuntimeError("No fertile mass specified!")

def setup_device(device):

    """ Run settings """
    device.settings.photon_transport = False
    device.settings.particles = int(1e3)
    device.settings.batches = 10

    # ==============================================================================
    # Tallies
    # ==============================================================================
    """ Cylindrical Mesh Tally """
    mesh = openmc.CylindricalMesh()
    mesh.r_grid = np.linspace(25, 200, num=25)
    mesh.z_grid = np.linspace(-200, 200, num=50)
    mesh.phi_grid = np.array([0, (2 * np.pi)/(18 * 2)])
    mesh_filter = openmc.MeshFilter(mesh)

    device.add_tally('Mesh Tally', ['flux', '(n,Xt)', 'heating-local', 'absorption'], filters=[mesh_filter])
    
    """ FLiBe Tally """
    device.add_tally('FLiBe Tally', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[])
    device.add_tally('Li Tally', ['(n,Xt)'], filters=[], nuclides=['Li6', 'Li7'])

    return device

# ==============================================================================
# Depletion Scan
# ==============================================================================

enrichments = [0, 1, 2, 4, 7.5, 12, 25, 50, 100]
np.savetxt(base_dir + '/enrichments.txt', enrichments)

for enrichment in enrichments:
    """ DEPLETION SETTINGS """
    print("~~~~~~~~~~~~~~~~~~ ENRICHMENT: " + str(enrichment) + " kg ~~~~~~~~~~~~~~~~~~")

    fusion_power = 500 #MW
    num_steps = 5
    time_steps = [100*24*60*60 / num_steps] * num_steps
    source_rates = [fusion_power * anp.neutrons_per_MJ] * num_steps

    chain_file = '/home/jlball/arc-nonproliferation/data/simple_chain_endfb71_pwr.xml'

    """ Generate blankets doped to specified enrichment """
    
    U_device = setup_device(anp.generate_device("U", mass, Li6_enrichment=enrichment))
    Th_device = setup_device(anp.generate_device("Th", mass, Li6_enrichment=enrichment))

    """ Run depletion calculation """

    # This is necessary to ensure the tally.xml file gets written into the directory where the depletion calculation will run
    os.mkdir(base_dir + '/Uranium/'+ str(enrichment))
    os.chdir(base_dir + '/Uranium/'+ str(enrichment))
    U_device.build()
    os.chdir('../../..')

    U_device.deplete(time_steps, 
        source_rates=source_rates, 
        operator_kwargs={'chain_file':chain_file, 'normalization_mode':'source-rate'}, 
        directory=base_dir + '/Uranium/'+ str(enrichment))

    os.mkdir(base_dir + '/Thorium/'+ str(enrichment))
    os.chdir(base_dir + '/Thorium/'+ str(enrichment))
    Th_device.build()
    os.chdir('../../..')

    Th_device.deplete(time_steps, 
        source_rates=source_rates, 
        operator_kwargs={'chain_file':chain_file, 'normalization_mode':'source-rate'}, 
        directory=base_dir + '/Thorium/' + str(enrichment))
