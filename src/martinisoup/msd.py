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

from .units import UnitSystem, unit_string


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
        Stop frame. Defaults to ``start + 100``, clipped to the trajectory length.
    fft : bool
        Use the FFT algorithm for MSD calculation. Default ``True``.

    Returns
    -------
    dict
        Keys:

        - ``resnames`` — list of unique metabolite type names
        - ``residue_timeseries`` — mean MSD per type (L²), shape ``(n_types, n_lagtimes)``
        - ``residue_std`` — standard deviation across particles (L²), shape ``(n_types, n_lagtimes)``
        - ``time`` — absolute times of the analysed frames (T)
        - ``lagtimes`` — lag time of each MSD point, ``index × dt`` (T)
        - ``msd`` — MSD averaged across all particles (L²)
        - ``dimensions`` — dimensionality factor used in the MSD calculation
        - ``dt`` — trajectory timestep (T)
        - ``mean_box_lengths`` — box lengths ``[a, b, c]`` averaged over the analysis frames (L)
        - ``units`` — pint unit string for each of the keys above

        L and T are the length and time units of the data as MDAnalysis holds
        it: Å and ps, as the Universe is read with default settings,
        regardless of the units used in the trajectory file. They are
        recorded explicitly in ``units``.
    """
    u = Universe(topology, trajectory)
    unit_system = UnitSystem.from_universe(u)

    n_frames = len(u.trajectory)
    start = start if start is not None else int(n_frames / 2)
    stop = stop if stop is not None else start + 100
    stop = min(stop, n_frames)

    # Average box lengths [a, b, c] (in the Universe's length unit) over the analysis
    # frames. Done before attaching NoJump so the transformation isn't applied
    # (and reset) needlessly while reading box dimensions.
    mean_box_lengths = np.mean(
        [ts.dimensions[:3] for ts in u.trajectory[start:stop]], axis=0
    )

    # NoJump must be attached on frame 0 so it initialises its reference correctly
    u.trajectory[0]
    u.trajectory.add_transformations(NoJump())

    MSD = msd_analysis.EinsteinMSD(u, select=metab_sel, fft=fft)
    MSD.run(start=start, stop=stop, verbose=True)

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
        'lagtimes': np.arange(len(msd_data['timeseries'])) * u.trajectory.dt,
        'msd': msd_data['timeseries'],
        'dimensions': MSD.dim_fac,
        'dt': u.trajectory.dt,
        'mean_box_lengths': mean_box_lengths,  # [a, b, c]
        'units': {
            'residue_timeseries': unit_string(unit_system.length_unit ** 2),
            'residue_std': unit_string(unit_system.length_unit ** 2),
            'time': unit_string(unit_system.time_unit),
            'lagtimes': unit_string(unit_system.time_unit),
            'msd': unit_string(unit_system.length_unit ** 2),
            'dimensions': 'dimensionless',
            'dt': unit_string(unit_system.time_unit),
            'mean_box_lengths': unit_string(unit_system.length_unit),
        },
    }
