# coding: utf8

import numpy as np
from ase.build import bulk

from msspec.calculator import MSSPEC, XRaySource
from msspec.iodata import Data
from msspec.utils import hemispherical_cluster, get_atom_index


# Create the cluster
a = 5.43
Si = bulk('Si', a=a, cubic=True)
cluster = hemispherical_cluster(Si,
                                diameter=30, planes=4,
                                emitter_plane=3,
                                shape = 'cylindrical',
                                )
for atom in cluster:
    atom.set('mean_square_vibration', 0.006)
    atom.set('mt_radius', 1.1)
cluster.emitter = get_atom_index(cluster, 0, 0, 0)

# Create a calculator and set parameters
calc = MSSPEC(spectroscopy='PED', algorithm='expansion')

calc.source_parameters.energy = XRaySource.AL_KALPHA
calc.source_parameters.theta  = -54.7
calc.source_parameters.phi    = 90
calc.spectroscopy_parameters.final_state = 1

calc.calculation_parameters.scattering_order = 3
calc.tmatrix_parameters.tl_threshold = 1e-4
calc.calculation_parameters.vibrational_damping = 'averaged_tl'
calc.calculation_parameters.RA_cutoff = 2
 
# Define path filtering options such that you only
# accept scattering paths with a forward cone <= 40°
# and whose length are <= cluster diameter
#
#

calc.set_atoms(cluster)

# Compute and add previous data for comparison
data = calc.get_theta_scan(level='2p',
                            theta=np.arange(-30., 80., 0.5),
                            phi=0,
                            kinetic_energy=1382.28)
no_filters = Data.load('path_filtering.hdf5')
data[0].add_columns(**{'no_filters': no_filters[0].cross_section})
view = data[0].views[0]
view.select('theta', 'cross_section', index=0, legend="With path filtering")
view.select('theta', 'no_filters', legend="Without path filtering")

data.view()
