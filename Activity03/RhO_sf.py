from ase import Atoms
from ase.visualize import view
from msspec.calculator import MSSPEC

# Create an atomic chain O-Rh
cluster = Atoms(['O', 'Rh'], positions = [(1,0,0), (0,0,4.)])

# Create a calculator
calc = MSSPEC(spectroscopy='PED')
calc.set_atoms(cluster)
cluster.emitter = 0

# Compute the scattering factor
data = calc.get_scattering_factors(kinetic_energy=723)

# Popup the results
data.view()
