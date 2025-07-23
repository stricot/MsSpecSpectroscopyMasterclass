from msspec.calculator import MSSPEC
from ase.build import fcc111, add_adsorbate
import numpy as np

data = None
all_z = ... # -> Define a list of z values for the adsorbate

for ... # -> Complete this for-loop over z values
    # construct the cluster
    cluster = fcc111('Rh', size = (2,2,1))
    cluster.pop(3)
    add_adsorbate(... # -> Put the oxygen atom on the fcc site
    cluster.emitter = ... # -> Oxygen is the last atom we added, so the indice is...

    # Define a calculator for single scattering calculations
    calc = MSSPEC(spectroscopy='PED', algorithm='expansion')
    calc.calculation_parameters.scattering_order = 1
    calc.set_atoms(cluster)

    # Compute
    data = calc.get_theta_phi_scan(level='1s', kinetic_energy=723, data=data)

data.view()
