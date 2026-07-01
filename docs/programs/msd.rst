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
     --stop INT            Stop frame (default: start + 100)
     --fft                 Use FFT algorithm for MSD calculation (default: True)
     --output PATH         Output pickle file (default: msd.pkl)

Output pickle keys:

- ``command`` — the full command used to produce the results
- ``resnames`` — list of unique metabolite type names
- ``residue_timeseries`` — mean MSD per type, shape ``(n_types, n_lagtimes)``
- ``residue_std`` — standard deviation across molecules, shape ``(n_types, n_lagtimes)``
- ``time`` — absolute times corresponding to each lagtime
- ``lagtimes`` — MSD values averaged across all molecules
- ``dimensions`` — dimensionality factor used in the MSD calculation
- ``dt`` — trajectory timestep in ps

----

msd-fitter
----------

Fits MSD output from ``martinisoup msd`` to a linear model and extracts diffusion
coefficients via the Einstein relation (D = slope / 6). Multiple input files are treated
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

Without ``--database``, outputs ``diffusion_coefficients.csv`` with columns ``resname``, ``D``, ``D_err``.
With ``--database``, results are grouped by class and the CSV has columns ``class``, ``resname``, ``D``, ``D_err``.

When multiple replica files are provided, the mean MSD across replicas is computed and the uncertainty
is propagated as ``sqrt(sum(σ_i²)) / n_replicas`` before fitting.
