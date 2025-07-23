from ase.build import mx2
from ase.visualize import view
from msspec.calculator import MSSPEC
from msspec.utils import hemispherical_cluster, get_atom_index

# Some usefull constants (a, c, d, D) for defining the structure
a=...

# Create the TiSe2 trilayer
# use ase help for this function
TiSe2 = mx2(formula=...)

# The preious cell is 2D, let's define the c-axis to take into account 
# the Van der Waals gap between trilayers
TiSe2.cell[2] = [0, 0, ...]

# To be aligned like in the paper
TiSe2.rotate(60, 'z', rotate_cell=True)

# Since the material is multi-elements, "tag" each inequivalent atom 
# of the unit cell with a number. The "Ti" atom is tagged 0 and "Se" 
# atoms are 1 and 2.
for i in range(3): 
    TiSe2[i].tag = i

cluster = hemispherical_cluster(TiSe2, emitter_tag=..., emitter_plane=..., planes=5)
cluster.emitter = get_atom_index(cluster, 0, 0, 0)

view(cluster)

# Create a calculator with Rehr-Albers series expansion algorithm
calc = MSSPEC(spectroscopy='PED', algorithm='expansion')
calc.set_atoms(cluster)

data = None
for ndif in range(1,4):    
    calc.calculation_parameters.scattering_order = ndif
    data = calc.get_theta_phi_scan(level='2p', kinetic_energy=1030, data=data)

data.view()