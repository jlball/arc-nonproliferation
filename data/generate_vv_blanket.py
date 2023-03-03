import numpy as np
import matplotlib.pyplot as plt
import sys

maj_raius = 400 #cm
min_radius = 120 #cm

elongation = 1.5
triangularity = 0.5 

blanket_thickness = 100 #cm

num_points = 20

def generate_point(theta, offset):
    radius = min_radius + offset
    x = maj_raius + radius*np.cos(theta + triangularity*np.sin(theta))
    y = elongation*radius*np.sin(theta)
    return np.array([x, y])

thetas = np.linspace(0, 2*np.pi, num=num_points)

vv_points = np.empty((num_points, 2))
blanket_points = np.empty((num_points, 2))
for i, theta in enumerate(thetas):
    vv_points[i] = generate_point(theta, 0)
    blanket_points[i] = generate_point(theta, blanket_thickness)

# Save points to text file
np.savetxt(str(sys.argv[1])+"_vv.txt", vv_points)
np.savetxt(str(sys.argv[1])+"_blanket.txt", blanket_points)

# Plot resulting shape
fig, ax = plt.subplots()
ax.plot(vv_points[:, 0], vv_points[:, 1])
ax.plot(blanket_points[:, 0], blanket_points[:, 1])
ax.set_aspect(1)
plt.show()