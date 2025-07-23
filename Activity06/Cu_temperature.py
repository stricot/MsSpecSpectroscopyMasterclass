from ase.build import bulk
import numpy as np

from msspec.calculator import MSSPEC, XRaySource
from msspec.utils import hemispherical_cluster, get_atom_index

def create_clusters(nplanes=3):
    copper  = bulk('Cu', a=3.6, cubic=True)
    clusters = []
    for emitter_plane in range(nplanes):
        cluster = hemispherical_cluster(copper,
                                        emitter_plane=emitter_plane,
                                        planes=emitter_plane+1,
                                        diameter=27,
                                        shape='cylindrical')
        cluster.absorber = get_atom_index(cluster, 0, 0, 0)
        # This is how to store extra information with your cluster
        cluster.info.update({
            'emitter_plane': emitter_plane,
        })
        clusters.append(cluster)
    return clusters


def compute(clusters, all_theta=[45., 83.],
            all_T=np.arange(300., 1000., 400.)):
    data = None
    for cluster in clusters:
        # Retrieve emitter's plane from cluster object
        plane   = cluster.info['emitter_plane']

        calc = MSSPEC(spectroscopy='PED', algorithm='expansion')
        calc.source_parameters.energy = XRaySource.AL_KALPHA
        calc.muffintin_parameters.interstitial_potential = 14.1

        # In simple scattering, it is common practice to use a real potential and
        # manually define a mean free path arbitrarily lower than the actual physical
        # value in an attempt to reproduce multiple scattering effects.
        calc.tmatrix_parameters.exchange_correlation = 'x_alpha_real'
        calc.calculation_parameters.mean_free_path = ... # -> half of the mean free
                                                         # path (see p6785)
        # Parameters for temperature effects
        calc.calculation_parameters.vibrational_damping = 'averaged_tl'
        calc.calculation_parameters.use_debye_model = ..... # Use the MsSpec help
        calc.calculation_parameters.debye_temperature = ... # and p6791 of the paper
        calc.calculation_parameters.vibration_scaling = ... # -> How much more do 
                                                            # surface atoms vibrate 
                                                            # than bulk atoms?
        calc.detector_parameters.average_sampling = 'low'
        calc.detector_parameters.angular_acceptance = 5.7

        calc.calculation_parameters.scattering_order = 1


        for T in all_T:
            # Define the sample temperature
            calc.calculation_parameters.temperature = T
            # Set the atoms and compute an azimuthal scan
            calc.set_atoms(cluster)
            data = calc.get_phi_scan(level='2p', theta=all_theta,
                                        phi=np.linspace(0, 100, 51),
                                        kinetic_energy=560, data=data)
            # Small changes to add some details in both the title of the dataset
            # and the figure
            view = data[-1].views[-1]
            t = view._plotopts['title'] + f" (plane #{plane:d}, T={T:.0f} K)"
            data[-1].title = t
            view.set_plot_options(autoscale=True, title=t)
        calc.shutdown()
    return data


def analysis(data, all_theta, all_T, nplanes):
    # Sum cross_section for all emitter's plane at a given T
    # Compute the anisotropy
    results = dict.fromkeys(all_T, []) 
    anisotropy = []
    for dset in data:
        # Retrieve temperature
        T = float(dset.get_parameter('CalculationParameters', 'temperature')['value'])
        # Update the sum in results
        if len(results[T]) == 0:
            results[T] = dset.cross_section
        else:
            results[T] += dset.cross_section

    anisotropy_dset = data.add_dset("Anisotropies")
    anisotropy_dset.add_columns(temperature=all_T)
    for theta in all_theta:
        col_name = f"theta{theta:.0f}"
        col_values = []
        i = np.where(dset.theta == theta)[0]
        for T in all_T:
            cs = results[T][i]
            Imax = np.max(cs)
            Imin = np.min(cs)
            A = (Imax - Imin)/Imax
            col_values.append(A)
        anisotropy_dset.add_columns(**{col_name:col_values/np.max(col_values)})

    
    anisotropy_view = anisotropy_dset.add_view('Anisotropies',
                         title='Relative anisotropies for Cu(2p)',
                         marker='o',
                         xlabel='T (K)',
                         ylabel=r'$\frac{\Delta I / I_{max}(T)}{\Delta I_{300}'
                                r'/ I_{max}(300)} (\%)$',
                         autoscale=True)
    for theta in all_theta:
        col_name = f"theta{theta:.0f}"
        anisotropy_view.select('temperature', col_name,
                                legend=r'$\theta = {:.0f} \degree$'.format(theta))

    return data



if __name__ == "__main__":
    nplanes = 4
    all_theta = np.array([45, 83])
    all_theta = np.array([300., 1000.])

    clusters = create_clusters(nplanes=nplanes)
    data = compute(clusters, all_T=all_T, all_theta=all_theta)
    data = analysis(data, all_T=all_T, all_theta=all_theta, nplanes=nplanes)
    data.view()
