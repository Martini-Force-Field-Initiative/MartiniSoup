msd / msd-fitter
=================

msd
---

Computes the mean squared displacement per metabolite type using the MDAnalysis Einstein
MSD implementation. A NoJump transformation is applied automatically to unwrap periodic
boundary conditions. Results are averaged across all copies of each molecule type.

.. code-block:: bash

   martinisoup msd topology.tpr trajectory.xtc --output msd.pkl

.. code-block:: text

   positional arguments:
     topology              Topology file (e.g. .tpr, .gro)
     trajectory            Trajectory file (e.g. .xtc, .trr)

   options:
     --metab_sel STR       MDAnalysis selection string for metabolites
     --start INT           Start frame (default: first frame of second half of trajectory)
     --stop INT            Stop frame (default: start + 100, clipped to trajectory length)
     --fft / --no-fft      Use FFT algorithm for MSD calculation (default: True)
     --output PATH         Output pickle file (default: msd.pkl)

Output pickle keys:

- ``command`` — the full command used to produce the results
- ``resnames`` — list of unique metabolite type names
- ``residue_timeseries`` — mean MSD per type (Å²), shape ``(n_types, n_lagtimes)``
- ``residue_std`` — standard deviation across particles (Å²), shape ``(n_types, n_lagtimes)``
- ``time`` — absolute times of the analysed frames (ps)
- ``lagtimes`` — lag time of each MSD point, ``index × dt`` (ps)
- ``msd`` — MSD averaged across all particles (Å²)
- ``dimensions`` — dimensionality factor used in the MSD calculation
- ``dt`` — trajectory timestep (ps)
- ``mean_box_lengths`` — box lengths ``[a, b, c]`` averaged over the analysis frames (Å)
- ``units`` — ``{key: pint unit string}`` for each of the keys above

Values are in the units of the data as MDAnalysis holds it, Å and ps, regardless of the
units used in the trajectory file (e.g. nm for ``.xtc``). The ``units`` entry records this
explicitly; see :doc:`../units`.

.. note::
   **Changed:** ``lagtimes`` now holds the lag times (ps). It previously held the MSD averaged across
   all particles (Å²), which is now under the new ``msd`` key. ``msd-fitter`` still reads
   files in the old layout and converts them automatically; other scripts that read
   ``lagtimes`` should switch to ``msd``. See :doc:`../changelog`.

----

msd-fitter
----------

Fits MSD output from ``martinisoup msd`` to a linear model and extracts diffusion
coefficients via the Einstein relation (D = slope / 2d, i.e. slope / 6 for 3D MSDs). Multiple input files are treated
as replicas and averaged before fitting.

.. code-block:: bash

   # per-molecule results (no database)
   martinisoup msd-fitter msd.pkl --output-dir results/

   # group by class using the remote M3-Metabolome database
   martinisoup msd-fitter msd.pkl --database --output-dir results/

   # group by class using a local database file
   martinisoup msd-fitter msd.pkl --database database.csv --output-dir results/

   # average over replicas before fitting
   martinisoup msd-fitter replica_1/msd.pkl replica_2/msd.pkl replica_3/msd.pkl --database --output-dir results/

.. code-block:: text

   positional arguments:
     files                 Pickle file(s) from `martinisoup msd`. Multiple files are averaged as replicas.

   options:
     --cut-start INT       Start index of the fitting window (default: 10)
     --cut-end INT         End index of the fitting window (default: 50)
     --database [PATH]     Group results by metabolite class. Use --database alone to fetch the
                           remote M3-Metabolome default, or supply a local CSV path.
                           Omit entirely for per-molecule results.
     --output-dir PATH     Directory for output plots and results CSV (default: current directory)
     --style PATH          Path to a matplotlib style file
     --temperature QUANTITY
                           Simulation temperature. When provided, applies the Yeh-Hummer
                           finite-size correction using the box dimensions stored in the MSD data.
                           A bare number is in K; a value with units, e.g. "37 degC", is converted.
     --eta QUANTITY        Dynamic viscosity for the Yeh-Hummer correction (default: 0.7103 mPa·s,
                           Martini 3 water). A bare number is in mPa·s; a value with units,
                           e.g. "0.00071 Pa*s", is converted.

Without ``--database``, outputs ``diffusion_coefficients.csv`` with columns ``resname``, ``D``, ``D_err``.
With ``--database``, results are grouped by class and the CSV has columns ``class``, ``resname``, ``D``, ``D_err``.
When ``--temperature`` is given, ``D`` is reported as ``D_PBC`` and an extra ``D_corr`` column holds the
Yeh-Hummer corrected value. All diffusion coefficients are in cm²/s.
``diffusion_coefficients.pkl`` holds ``command``, ``results`` and
``units`` (``{'results': 'centimeter ** 2 / second'}``).

Files in the old layout (MSD values under ``lagtimes``, no ``msd`` key) are recognised and
converted, with lag times rebuilt from ``dt``.

Each input file is converted from its recorded units to Å and ps before fitting, so replicas
saved in different units are combined correctly. Files from older versions without recorded
units are assumed to be in Å and ps (with a warning); their ``trajectory_units`` entry, which
describes the trajectory file rather than the data, is ignored.

When multiple replica files are provided, the mean MSD across replicas is computed and the uncertainty
is propagated as ``sqrt(sum(σ_i²)) / n_replicas`` before fitting.

----

Method notes
------------

Two aspects of the current method are worth bearing in mind when interpreting results.

**MSDs are computed per bead, not per molecule.**
``martinisoup msd`` runs the MDAnalysis Einstein MSD on every bead in the selection and then averages
over all beads of a given residue name. It does not use molecular centres of mass. At long lag times
the bead and centre-of-mass MSDs have the same slope, so ``D`` is unaffected. At short lag times,
however, internal flexibility and rotation of multi-bead molecules add to the bead MSD, raising the
curve above the centre-of-mass MSD. This is another reason to keep the fitting window away from the
shortest lag times. In addition, ``residue_std`` is the spread across beads, which are correlated
within a molecule, rather than across independent molecules. The resulting uncertainties (and the
fit weights derived from them) are therefore optimistic, particularly for molecules with many beads.

**The fitting window is set in frames, not time.**
``--cut-start`` and ``--cut-end`` are frame indices into the lag-time axis, so the time range they
cover is ``index × dt`` and depends on how often the trajectory was written. The same defaults
(10–50) cover 100–500 ps for a trajectory saved every 10 ps, but 1–5 ns for one saved every 100 ps.
If the window falls in the short-time (ballistic or subdiffusive) regime, the slope, and hence ``D``,
will be underestimated. At the other end, MSD points near the end of the analysed range are averaged
over few time origins and are noisy. Always check the per-molecule plots to confirm that the fitted
window lies in the linear part of the curve, and adjust the indices for your output frequency.
