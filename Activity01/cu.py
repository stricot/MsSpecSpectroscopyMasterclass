from ase.io import read
from msspec.calculator import MSSPEC


cluster = read('copper.cif')
# view the cluster
cluster.edit()

# The "emitter" atom is located in the middle of the 3rd plane
cluster.emitter = 10

# Create a "calculator"
calc = MSSPEC(spectroscopy='PED', algorithm='inversion')
calc.set_atoms(cluster)

data = calc.get_theta_scan(level='2p3/2')

# Plot the result with the interactive GUI
data.view()

# Or plot using matplotlib directly
from matplotlib import pyplot as plt
data[0].views[0].plot()
plt.show()
