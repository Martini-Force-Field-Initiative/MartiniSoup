residence-times / residence-fitter
====================================

residence-times
----------------

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

residence-fitter
------------------

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
