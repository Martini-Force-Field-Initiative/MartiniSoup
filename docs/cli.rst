CLI Reference
=============

All commands are available as subcommands of the ``martinisoup`` CLI:

.. code-block:: text

   martinisoup <command> [options]

   Commands:
     binding-frequency    Metabolite–protein contact frequency analysis
     residence-times      Binding event residence time analysis
     msd                  Mean squared displacement per metabolite type
     msd-fitter           Fit MSD data to extract diffusion coefficients
     residence-fitter     Fit residence-time histograms by metabolite class

Default selections
------------------

Unless overridden, MDAnalysis selection strings default to:

.. code-block:: text

   metabolite_selection = 'not resname NA CL ION GLU ASN VAL ALA GLY ARG SER PRO THR PHE GLN LYS LEU ASP ILE MET HIS CYS TRP TYR'
   protein_selection    = 'resname GLU ASN VAL ALA GLY ARG SER PRO THR PHE GLN LYS LEU ASP ILE MET HIS CYS TRP TYR'

----

binding-frequency
-----------------

Counts how often each metabolite type contacts each protein type over a trajectory.
Results are normalised per protein monomer per metabolite copy.

.. code-block:: bash

   martinisoup binding-frequency topology.tpr trajectory.xtc --output results.pkl

.. code-block:: text

   positional arguments:
     topology              Topology file (e.g. .pdb, .tpr)
     trajectory            Trajectory file (e.g. .xtc)

   options:
     --protein_sel STR     MDAnalysis selection string for proteins
     --metab_sel STR       MDAnalysis selection string for metabolites
     --cutoff FLOAT        Distance cutoff in Å (default: 5.0)
     --parallel            Run in parallel mode
     --n_workers INT       Number of parallel workers (default: 4)
     --chunk_size INT      Frames per parallel chunk (default: 100)
     --output PATH         Path to save results as a pickle file (default: binding.pkl)

Output pickle keys:

- ``command`` — the full command used to produce the results
- ``protein_results`` — per-protein, per-metabolite raw and normalised contact counts
- ``n_metabolites`` — metabolite copy numbers in the system
- ``n_proteins`` — protein monomer counts in the system

----

residence-times
---------------

Tracks binding and unbinding events for each metabolite molecule across the trajectory
and records event durations.

.. note::
   The ``--summary`` flag is strongly recommended. Raw output can be very large.

.. code-block:: bash

   martinisoup residence-times topology.tpr trajectory.xtc --output lifetimes.pkl --summary summary.pkl

.. code-block:: text

   positional arguments:
     topology              Topology file (e.g. .pdb, .gro, .tpr)
     trajectory            Trajectory file (e.g. .xtc, .dcd, .trr)

   options:
     --metabolite-selection STR   MDAnalysis selection string for metabolites
     --protein-selection STR      MDAnalysis selection string for proteins
     --cutoff FLOAT               Distance cutoff in Å (default: 5.0)
     --start INT                  Start frame (default: 0)
     --stop INT                   Stop frame (default: end of trajectory)
     --step INT                   Frame step size (default: 1)
     --output PATH                Output pickle file (default: lifetimes.pkl)
     --summary FILE               Write a summary pickle (unique residence times and
                                  their counts) to FILE
     --parallel                   Run contact detection in parallel (Phase 1), then
                                  run the state machine sequentially (Phase 2).
                                  Produces identical results to serial mode.
     --n_workers INT              Number of worker processes for parallel mode (default: 4)
     --chunk_size INT             Frames per chunk in parallel mode (default: 100)

Output pickle keys:

- ``command`` — the full command used to produce the results
- ``residences`` — ``{molecule_type: [durations]}`` dictionary

----

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

----

residence-fitter
----------------

Fits residence-time data from ``martinisoup residence-times`` to a power law model.
Multiple input files are treated as replicas and averaged before fitting.

.. code-block:: bash

   # per-molecule results (no database)
   martinisoup residence-fitter lifetimes_summary.pkl --output-dir results/

   # group by class using the remote M3-Metabolome database
   martinisoup residence-fitter lifetimes_summary.pkl --database --output-dir results/

   # group by class using a local database file and custom fit weights
   martinisoup residence-fitter lifetimes_summary.pkl --database database.csv --weights weights.json --output-dir results/

   # raw (unsummarised) output, averaged over replicas
   martinisoup residence-fitter replica_1/lifetimes.pkl replica_2/lifetimes.pkl --unsummarised --database --output-dir results/

.. code-block:: text

   positional arguments:
     files                 Pickle file(s) from `martinisoup residence-times`. Multiple files are averaged as replicas.

   options:
     --database [PATH]     Group results by metabolite class. Use --database alone to fetch the
                           remote M3-Metabolome default, or supply a local CSV path.
                           Omit entirely for per-molecule results.
     --unsummarised        Input files are raw (unsummarised) residence-times output. By default,
                           summarised output (produced with --summary) is expected.
     --bins-start FLOAT    Start of log-spaced bin range as log10 value in ns (default: 0)
     --bins-stop FLOAT     End of log-spaced bin range as log10 value in ns (default: log10(500))
     --bins-n INT          Number of bins (default: 25)
     --weights PATH        JSON file mapping class names to fit weight upper bounds
     --output-dir PATH     Directory for output plots and results CSV (default: current directory)
     --style PATH          Path to a matplotlib style file

Outputs:

- ``residence_fits.png`` — grid figure with one subplot per metabolite class showing the power law fit
- ``residence_exponents.png`` — bar chart of fitted exponents
- ``residence_exponents.csv`` — columns ``class``, ``exponent``, ``exponent_err``

The ``--weights`` JSON maps class names to the upper bound of the log-spaced fit weights:

.. code-block:: json

   {
     "Ions": 500,
     "Nucleotides": 10000
   }
