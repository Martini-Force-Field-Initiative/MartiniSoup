"""
Mean squared displacement analysis per metabolite type.

Uses the MDAnalysis Einstein MSD implementation with a NoJump transformation
to handle periodic boundary conditions. Results are averaged across all copies
of each molecule type.
"""

from __future__ import annotations

import numpy as np
from tqdm import tqdm
from MDAnalysis import Universe
from MDAnalysis.analysis import msd as msd_analysis
from MDAnalysis.transformations import NoJump


def compute_msd(
    topology: str,
    trajectory: str,
    metab_sel: str,
    start: int | None = None,
    stop: int | None = None,
    fft: bool = True,
) -> dict:
    """
    Compute the mean squared displacement per metabolite type.

    A NoJump transformation is applied automatically to unwrap periodic
    boundary conditions. Results are averaged across all copies of each
    molecule type.

    Parameters
    ----------
    topology : str
        Path to a topology file (e.g. ``.tpr``, ``.gro``).
    trajectory : str
        Path to a trajectory file (e.g. ``.xtc``, ``.trr``).
    metab_sel : str
        MDAnalysis selection string for metabolites.
    start : int, optional
        Start frame. Defaults to the first frame of the second half of the
        trajectory.
    stop : int, optional
        Stop frame. Defaults to ``start + 100``.
    fft : bool
        Use the FFT algorithm for MSD calculation. Default ``True``.

    Returns
    -------
    dict
        Keys:

        - ``resnames`` — list of unique metabolite type names
        - ``residue_timeseries`` — mean MSD per type, shape ``(n_types, n_lagtimes)``
        - ``residue_std`` — standard deviation across molecules, shape ``(n_types, n_lagtimes)``
        - ``time`` — absolute times corresponding to each lagtime
        - ``lagtimes`` — MSD values averaged across all molecules
        - ``dimensions`` — dimensionality factor used in the MSD calculation
        - ``dt`` — trajectory timestep in ps
    """
    u = Universe(topology, trajectory)
    u.trajectory.add_transformations(NoJump())

    start = start if start is not None else int(len(u.trajectory) / 2)
    stop = stop if stop is not None else start + 100

    MSD = msd_analysis.EinsteinMSD(u, select=metab_sel, fft=fft)
    MSD.run(start=start, stop=stop, verbose=True)

    # Average box lengths [a, b, c] in Å (MDAnalysis internal) over the analysis frames
    mean_box_lengths = np.mean(
        [u.trajectory[i].dimensions[:3] for i in range(start, stop)], axis=0
    )

    msd_data = dict(MSD.results)
    msd_array = msd_data['msds_by_particle']

    metabolites = u.select_atoms(metab_sel)
    unique_resnames = list(set(metabolites.resnames))

    residue_average = np.zeros((len(unique_resnames), msd_array.shape[0]))
    residue_std = np.zeros((len(unique_resnames), msd_array.shape[0]))

    for i, resname in enumerate(tqdm(unique_resnames, desc="Averaging per molecule type")):
        mask = metabolites.resnames == resname
        residue_average[i] = np.mean(msd_array.T[mask], axis=0)
        residue_std[i] = np.std(msd_array.T[mask], axis=0)

    return {
        'resnames': unique_resnames,
        'residue_timeseries': residue_average,
        'residue_std': residue_std,
        'time': MSD.times,
        'lagtimes': msd_data['timeseries'],
        'dimensions': MSD.dim_fac,
        'dt': u.trajectory.dt,
        'trajectory_units': dict(u.trajectory.units),
        'mean_box_lengths': mean_box_lengths,  # [a, b, c] in Å
    }
