import arc_nonproliferation as anp
import openmc
import openmc.deplete
import numpy as np
import os
import sys
import pickle

""" Handle command line arguments and setup directory structure """
if sys.argv[1] is not None:
    base_dir = str(sys.argv[1])
    os.mkdir(base_dir)
    print("OpenMC output saved to:", base_dir)

    os.mkdir(base_dir + '/Uranium')
    os.mkdir(base_dir + '/Thorium')

# This function handles the simulation specific setup of each device object
def setup_device(device):
    """ Run settings """
    device.settings.photon_transport = False
    device.settings.particles = int(1e5)
    device.settings.batches = 10
    device.settings.survival_biasing = True

    """ Cylindrical Mesh Tally """
    mesh = openmc.CylindricalMesh()
    mesh.r_grid = np.linspace(100, 700, num=50)
    mesh.z_grid = np.linspace(-300, 300, num=50)
    mesh.phi_grid = np.array([0, (2 * np.pi)/(18 * 2)])
    mesh_filter = openmc.MeshFilter(mesh)

    # """ Spectrum Mesh """
    # spec_mesh = openmc.CylindricalMesh()
    # spec_mesh.r_grid = np.linspace(100, 700, num=50)
    # spec_mesh.z_grid = np.linspace(-15, 15, num=1)
    # spec_mesh.phi_grid = np.array([0, (2 * np.pi)/(18 * 2)])
    # spec_mesh_filter = openmc.MeshFilter(mesh)

    """ Cell Filter """
    blanket_cell = device.get_cell(name='blanket')
    blanket_filter = openmc.CellFilter(blanket_cell)

    """ Energy Filter """
    energy_filter = openmc.EnergyFilter.from_group_structure("CCFE-709")

    device.add_tally('Mesh Tally', ['flux', '(n,Xt)', 'heating-local', 'absorption'], filters=[mesh_filter])

    """ FLiBe Tally """
    #flibe_filter = openmc.MaterialFilter(anp.get_material_by_name(device.materials, "doped_flibe"))
    device.add_tally('FLiBe Tally', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[])
    device.add_tally('Flux Tally', ['flux'], filters=[energy_filter, blanket_filter])
    device.add_tally('Li Tally', ['(n,Xt)'], filters=[], nuclides=['Li6', 'Li7'])
    device.add_tally("Feritle Tally", ['absorption', '(n,gamma)', 'fission'], nuclides=['U238, U235'], filters=[])

    return device

# ==============================================================================
# Depletion Run
# ==============================================================================

masses = np.array([5e3, 7e3, 10e3, 15e3, 30e3, 50e3]) #kg of fertile material
np.savetxt(base_dir + '/masses.txt', masses)

""" DEPLETION SETTINGS """
for mass in masses:
    print("~~~~~~~~~~~~~~~~~~ FERTILE MASS: " + str(mass) + " kg ~~~~~~~~~~~~~~~~~~")

    fusion_power = 500 #MW
    num_steps = 50
    time_steps = [365 / num_steps] * num_steps
    source_rates = [fusion_power * anp.neutrons_per_MJ] * num_steps

    chain_file = '/home/jlball/arc-nonproliferation/data/chain_endfb71_pwr.xml'
    openmc.config['chain_file'] = chain_file

    """ Generate blankets doped to specified mass """

    U_device = setup_device(anp.generate_device("U", mass))
    Th_device = setup_device(anp.generate_device("Th", mass))

    """ Run depletion calculation """
    # This is necessary to ensure the tally.xml file gets written into the directory where the depletion calculation will run
    os.mkdir(base_dir + '/Uranium/'+ str(mass))
    os.chdir(base_dir + '/Uranium/'+ str(mass))
    U_device.build()

    U_flux, U_micro_xs = openmc.deplete.get_microxs_and_flux(U_device,
                                                            [U_device.channel, U_device.blanket],
                                                            run_kwargs = {"threads":20,
                                                                        "particles":int(1e3)})

    flux_file = open('U_flux', 'ab')
    pickle.dump(U_flux, flux_file)

    microxs_file = open("U_microxs", 'ab')
    pickle.dump(U_micro_xs, microxs_file)

    U_operator = openmc.deplete.IndependentOperator(openmc.Materials([U_device.channel.fill, U_device.blanket.fill]), 
                                                    U_flux,
                                                    U_micro_xs,
                                                    chain_file=chain_file, 
                                                    normalization_mode='source-rate', 
                                                    reduce_chain=True, 
                                                    reduce_chain_level=5, 
                                                    )

    U_integrator = openmc.deplete.PredictorIntegrator(U_operator, 
                                                    time_steps, 
                                                    source_rates=source_rates,
                                                    timestep_units='d')

    U_integrator.integrate()

    os.chdir('../../..')

    os.mkdir(base_dir + '/Thorium/'+ str(mass))
    os.chdir(base_dir + '/Thorium/'+ str(mass))
    Th_device.build()

    Th_flux, Th_micro_xs = openmc.deplete.get_microxs_and_flux(Th_device,
                                                                    [Th_device.channel, Th_device.blanket],
                                                                    run_kwargs = {"threads":20,
                                                                                    "particles":int(1e3)})

    Th_operator = openmc.deplete.IndependentOperator(openmc.Materials([Th_device.channel.fill, Th_device.blanket.fill]), 
                                                    Th_flux,
                                                    Th_micro_xs,
                                                    chain_file=chain_file, 
                                                    normalization_mode='source-rate', 
                                                    reduce_chain=True, 
                                                    reduce_chain_level=5, 
                                                    )

    Th_integrator = openmc.deplete.PredictorIntegrator(Th_operator, 
                                                    time_steps, 
                                                    source_rates=source_rates,
                                                    timestep_units='d')

    Th_integrator.integrate()

    os.chdir('../../..')

