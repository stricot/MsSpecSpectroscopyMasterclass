from ase import Atoms
from ase.build import add_adsorbate, bulk

from msspec.calculator import MSSPEC, RFACTOR
from msspec.utils import hemispherical_cluster
from msspec.looper import Sweep, Looper

import numpy as np


def create_cluster(height=1., theta=45, phi=0, bond_length=1.15):
    # Fill the body of this function. The 'cluster' object in built according to
    # values provided by the keyword arguments:
    # height (in angströms): the 'z' distance between the Fe surface and the C atom
    # theta and phi (in degrees): the polar and azimuthal orientation of the CP molecule
    #                             (theta=0° aligns the molecule withe the surface normal
    #                              phi=0° corresponds to the [100] direction of iron)
    # bond_length (in angströms): the C-O distance

    iron = bulk('Fe', cubic=True)
    cluster = hemispherical_cluster(iron, diameter=5, planes=2, emitter_plane=1)

    t = np.radians(theta)
    p = np.radians(phi)

    z = bond_length * np.cos(t)
    x = bond_length * np.sin(t) * np.cos(p)
    y = bond_length * np.sin(t) * np.sin(p)
    CO=Atoms('CO',positions=[(0,0,0),(x,y,z)])

    add_adsorbate(cluster,CO, height=height)

    # Keep those 2 lines at the end of your function
    # Store some information in the cluster object
    cluster.info.update(adsorbate={'theta': theta, 'phi': phi, 'height': height, 
                                   'bond_length': bond_length})
    return cluster


def compute_polar_scan(cluster, folder='calc'):
    calc = MSSPEC(spectroscopy='PED', algorithm='expansion', folder=folder)
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
    calc.shutdown()

    return data


###############################################################################
# Main part
###############################################################################
# 1) Multiprocess calculations 
theta = Sweep(key='theta', comments="The molecule tilt angle",
              start=50, stop=60, step=1, unit='degree')
phi = Sweep(key='phi', comments="The molecule azimuthal angle",
            values=[0,45], unit='degree')

def process(theta, phi, **kwargs):
    cluster = create_cluster(theta=theta, phi=phi, height=0.6, bond_length=1.157)
    i = kwargs.get('sweep_index')
    data = compute_polar_scan(cluster, folder=f'calc_{i:d}')
    dset = data[-1]
    return dset.theta, dset.cross_section

looper = Looper()
looper.pipeline = process
df = looper.run(theta, phi, ncpu=4)

# Black magic to convert the pandas dataframe object 'df' to the 
# parameters dict and the resulst list (will be easier in a future
# version ;-) ).
parameters = df.to_dict('list')
results = np.reshape(parameters.pop('output'), (df.shape[0]*2,-1))

# 2) R-Factor analysis
# Load the experimental data
exp_data = np.loadtxt('experimental_data.txt')

# Create an R-Factor calculator
rfc = RFACTOR()
rfc.set_references(exp_data[:,0], exp_data[:,1])

# Perform the R-Factor analysis
data = rfc.run(*results, **parameters)
data.view()
