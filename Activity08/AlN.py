from ase.build import bulk
import numpy as np
from msspec.calculator import MSSPEC, XRaySource
from msspec.utils import hemispherical_cluster, get_atom_index

def create_clusters(nplanes=6):
    def get_AlN_tags_planes(side, emitter):
        AlN = # AlN is a Wurtzite crystal with a=3.11 and c=4.975 angstroms       # <= HERE
        [atom.set('tag', i) for i, atom in enumerate(AlN)]
        if side == 'Al':
            AlN.rotate([0,0,1],[0,0,-1])
            Al_planes = range(0, nplanes, 2)
            N_planes  = range(1, nplanes, 2)
        else:
            N_planes  = range(0, nplanes, 2)
            Al_planes = range(1, nplanes, 2)
        if emitter == 'Al':
            tags = [0, 2]
            planes = Al_planes
        else:
           tags = [1, 3]
           planes = N_planes
        return AlN, tags, planes

    clusters = []
    for side in ('Al', 'N'):
        for emitter in ('Al', 'N'):
            AlN, tags, planes = get_AlN_tags_planes(side, emitter)
            for emitter_tag in tags:
                for emitter_plane in planes:
                    cluster = # hemis…, construct the cluster here with           # <= HERE
                              # 2 planes below the emitter
                    cluster.absorber = get_atom_index(cluster, 0, 0, 0)
                    cluster.info.update({
                        'emitter_plane': emitter_plane,
                        'emitter_tag'  : emitter_tag,
                        'emitter'      : emitter,
                        'side'         : side,
                    })
                    clusters.append(cluster)
                    print("Added cluster {}-side, emitter {}(tag {:d}) in "
                          "plane #{:d}".format(side, emitter, emitter_tag,
                                               emitter_plane))
    return clusters


def compute(clusters, theta=np.arange(-20., 80., 1.), phi=0.):
    data = None
    for ic, cluster in enumerate(clusters):
        # Retrieve info from cluster object
        side    = cluster.info['side']
        emitter = cluster.info['emitter']
        plane   = cluster.info['emitter_plane']
        tag     = cluster.info['emitter_tag']

        # Set the level and the kinetic energy
        if emitter == 'Al':
            level = #####                                                         # <= HERE
            ke    = #####                                                         # <= HERE
        elif emitter == 'N':
            level = #####                                                         # <= HERE
            ke    = #####                                                         # <= HERE

        calc = # Create a calculator using the RA series expansion algorithm      # <= HERE

        calc.source_parameters.energy = #####                                     # <= HERE
        calc.source_parameters.theta  = #####                                     # <= HERE

        calc.detector_parameters.angular_acceptance = #####                       # <= HERE
        calc.detector_parameters.average_sampling   = 'medium'

        calc.calculation_parameters.scattering_order = max(1, min(4, plane))
        calc.calculation_parameters.path_filtering  = 'forward_scattering'
        calc.calculation_parameters.off_cone_events = 1
        [a.set('forward_angle', 30.) for a in cluster]

        calc.set_atoms(cluster)

        data = calc.get_theta_scan(level=level, theta=theta, phi=phi,
                                   kinetic_energy=ke, data=data)
        dset = data[-1]
        dset.title = "\'{}\' side - {}({}) tag #{:d}, plane #{:d}".format(
            side, emitter, level, tag, plane)

    return data


def analysis(data):
    tmp_data = {}
    for dset in data:
        info = dset.get_cluster().info
        side = info['side']
        emitter = info['emitter']
        try:
            key = '{}_{}'.format(side, emitter)
            tmp_data[key] += dset.cross_section
        except KeyError:
            tmp_data[key] = dset.cross_section.copy()

    tmp_data['theta']   = dset.theta.copy()
    tmp_data['Al_side'] = tmp_data['Al_Al'] / tmp_data['Al_N']
    tmp_data['N_side']  = tmp_data['N_Al']  / tmp_data['N_N']

    # now add all columns
    substrate_dset = data.add_dset('Total substrate signal')
    substrate_dset.add_columns(**tmp_data)

    view = substrate_dset.add_view('Ratios',
                                   title=r'Al(2p)/N(1s) ratios on both polar '
                                         r'sides of AlN in the (10$\bar{1}$0) '
                                         r'azimuthal plane',
                                   xlabel=r'$\Theta (\degree$)',
                                   ylabel='Intenisty ratio')
    view.select('theta', 'Al_side', legend='Al side',
                where="theta >= 0 and theta <=70")
    view.select('theta', 'N_side', legend='N side',
                where="theta >= 0 and theta <=70")
    view.set_plot_options(autoscale=True)

    return data


clusters = create_clusters()
data     = compute(clusters)
data     = analysis(data)
data.view()

