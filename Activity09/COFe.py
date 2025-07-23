from ase import Atoms
from ase.build import add_adsorbate, bulk

from msspec.calculator import MSSPEC, RFACTOR
from msspec.utils import hemispherical_cluster

import numpy as np


def create_cluster(height=1., theta=45, phi=0, bond_length=1.15):
    # Fill the body of this function. The 'cluster' object in built according to
    # values provided by the keyword arguments:
    # height (in angströms): the 'z' distance between the Fe surface and the C atom
    # theta and phi (in degrees): the polar and azimuthal orientation of the CP molecule
    #                             (theta=0° aligns the molecule withe the surface normal
    #                              phi=0° corresponds to the [100] direction of iron)
    # bond_length (in angströms): the C-O distance

    # Keep those 2 lines at the end of your function
    # Store some information in the cluster object
    cluster.info.update(adsorbate={'theta': theta, 'phi': phi, 'height': height, 'bond_length': bond_length})
    return cluster


def compute_polar_scan(cluster):
    calc = MSSPEC(spectroscopy='PED', algorithm='expansion')
    calc.set_atoms(cluster)

    # SSC calculations
    calc.calculation_parameters.scattering_order = 1

    # Add temperature effects
    [atom.set('mean_square_vibration', 0.005) for atom in cluster]
    calc.calculation_parameters.vibrational_damping = 'averaged_tl'

    polar_angles = np.arange(-5, 85, 0.5)
    # set the Carbon as absorber and compute the polar scan
    cluster.absorber = cluster.get_chemical_symbols().index('C')
    data = calc.get_theta_scan(level='1s', theta=polar_angles, kinetic_energy=1202)

    return data


###############################################################################
# Main part
###############################################################################
results    = [] # polar angles and calculated cross_sections will be appended
                # to this list after each 'compute_polar_scan' call
parameters = {'theta': [], 'phi': []} # and corresponding parameters will also
                                      # be stored in this dictionary

# 1) Run calculations for different geometries
for theta in ...
    for phi in ...
        # Create the cluster
        cluster = ...

        # Compute
        data = ...

        # Update lists of results and parameters
        results.append(data[-1].theta.copy())
        results.append(data[-1].cross_section.copy())
        parameters['theta'].append(theta)
        parameters['phi'].append(phi)

# 2) R-Factor analysis
# Load the experimental data
exp_data = np.loadtxt('experimental_data.txt')

# Create an R-Factor calculator
rfc = RFACTOR()
rfc.set_references(exp_data[:,0], exp_data[:,1])

# Perform the R-Factor analysis
data = rfc.run(*results, **parameters)
data.view()
