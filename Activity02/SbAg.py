from ase.build import bulk
from ase.visualize import view

from msspec.calculator import MSSPEC
from msspec.utils import hemispherical_cluster, get_atom_index, cut_plane
import numpy as np
from matplotlib import pyplot as plt

# Create the silver cell
Ag = bulk('Ag', cubic=True)
# Orientate the cell in the [111] direction
Ag.rotate((1,1,1), (0,0,1), rotate_cell=True)
# Align the azimuth to match experimental reference
Ag.rotate(15, 'z', rotate_cell=True)

# Create a cluster
cluster = hemispherical_cluster(Ag, diameter=20, emitter_plane=0)
cluster = cut_plane(cluster, z=-4.8)
cluster.emitter = get_atom_index(cluster, 0,0,0)
cluster[cluster.emitter].symbol = 'Sb'

# Create a calculator
calc = MSSPEC(spectroscopy='PED', algorithm='inversion')
calc.set_atoms(cluster)

# Define parameters
calc.source_parameters.theta = 0
calc.source_parameters.phi = 0

calc.detector_parameters.angular_acceptance = 1
calc.detector_parameters.average_sampling = 'low'

calc.muffintin_parameters.interstitial_potential = 0

# Compute an azimuthal scan
data = calc.get_phi_scan(level='4d', theta=40, phi=np.linspace(0,240,121), kinetic_energy=45)

# Normalize data between [0,1] (to ease comparison with experimental data)
dset = data[0]
dset.cross_section -= dset.cross_section.min()
dset.cross_section /= dset.cross_section.max()

# Add experimental data points in the dataset
x, y = np.loadtxt('data.txt').T
dset.add_columns(experiment=y)

# Add points to view
view = dset.views[0]
view.select('phi', 'experiment', legend='Exp. data')

# Popup GUI
data.view()

# Remove temp. files
calc.shutdown()
