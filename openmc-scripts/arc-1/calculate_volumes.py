import openmc
import arc_nonproliferation as anp

# Generate arbitrary device
device = anp.generate_device('U', 0)

device.settings.photon_transport = False
device.settings.particles = int(1e3)
device.settings.batches = 5

lower_left = (-650, -650, -350)
upper_right = (650, 650, 350)

blanket = device.get_cell('blanket')
channels = device.get_cell('channels')

vol_calc = openmc.VolumeCalculation([blanket, channels], int(1e6),
                                    lower_left, upper_right)

device.settings.volume_calculations = [vol_calc]

vol_calc.set_trigger(1e-02, 'rel_err')

device.build()
openmc.calculate_volumes()

#vol_calc_load = openmc.VolumeCalculation.from_hdf5('volume_1.h5')
#print(vol_calc_load.volumes[7].n)