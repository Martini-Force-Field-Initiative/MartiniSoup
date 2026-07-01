Python API
==========

The public API is importable directly from the ``martinisoup`` package:

.. code-block:: python

   from martinisoup import (
       track_serial, track_parallel, compute_msd, BindingState, ResidenceRegistry,
       compute_rdf, compute_rdf_parallel,
   )

Parallel execution
-------------------

Every ``*_parallel`` function in the package (``compute_rdf_parallel``,
``analyse_trajectory_parallel`` in both ``protein_clustering`` and
``metabolite_partitioning``, ``count_contacts_parallel``, and the geometry
pass of ``track_parallel``) is built on a single shared engine in
``martinisoup.parallel``. The trajectory is split into frame-index chunks,
each chunk is dispatched to a worker process that rebuilds the
``MDAnalysis.Universe`` from the topology/trajectory file paths (Universes
cannot be pickled across processes) and steps through its assigned frames,
and results are reassembled in ascending frame order in the main process.

A given analysis only needs to supply two module-level functions:

- ``setup_fn(u, *setup_args) -> context`` — run once per worker after the
  Universe is rebuilt (e.g. atom selections, static per-topology metadata).
- ``per_frame_fn(u, context, ts) -> result`` — run once per frame.

Both must be importable module-level functions rather than closures or
lambdas, since they are pickled to worker processes as part of the task
arguments.

.. autofunction:: martinisoup.parallel.map_trajectory_parallel

.. autofunction:: martinisoup.parallel.chunk_frame_indices

MSD analysis
------------

.. autofunction:: martinisoup.compute_msd

MSD fitting
-----------

These functions underlie the ``msd-fitter`` command.

.. autofunction:: martinisoup.scripts.msd_fitter.load_datasets

.. autofunction:: martinisoup.scripts.msd_fitter.build_lagtimes

.. autofunction:: martinisoup.scripts.msd_fitter.average_replicas

.. autofunction:: martinisoup.scripts.msd_fitter.fit_and_plot

.. autofunction:: martinisoup.scripts.msd_fitter.save_results

Residence time tracking
-----------------------

``track_parallel`` only parallelises the per-frame contact query (finding
metabolite atoms within cutoff of a protein atom); the binding/unbinding
state machine that turns those per-frame contacts into durations is
inherently sequential (a binding event can span many frames) and runs once,
in the main process, via the same fold used by ``track_serial``.

.. autofunction:: martinisoup.track_serial

.. autofunction:: martinisoup.track_parallel

Residence time fitting
----------------------

These functions underlie the ``residence-fitter`` command.

.. autofunction:: martinisoup.scripts.residence_fitter.load_datasets

.. autofunction:: martinisoup.scripts.residence_fitter.build_histograms

.. autofunction:: martinisoup.scripts.residence_fitter.fit_and_plot

.. autofunction:: martinisoup.scripts.residence_fitter.save_results

Data structures
---------------

.. autoclass:: martinisoup.BindingState
   :members:

.. autoclass:: martinisoup.ResidenceRegistry
   :members:
   :undoc-members:

Contact analysis
----------------

These functions underlie the ``binding-frequency`` command and can be called
directly for custom workflows. ``count_contacts_parallel`` accepts the same
``start``/``stop``/``step`` frame-range parameters as ``count_contacts``
(previously it always processed the whole trajectory); the CLI does not yet
expose these for ``binding-frequency``.

.. autofunction:: martinisoup.contact.count_contacts

.. autofunction:: martinisoup.contact.count_contacts_parallel

.. autofunction:: martinisoup.scripts.binding_frequency.index_protein_segments

Database
--------

.. autofunction:: martinisoup.database.load_metabolite_classes

Metabolite partitioning
------------------------

These functions underlie the ``metabolite-partitioning`` command.

.. autofunction:: martinisoup.metabolite_partitioning.analyse_frame

.. autofunction:: martinisoup.metabolite_partitioning.analyse_trajectory

.. autofunction:: martinisoup.metabolite_partitioning.analyse_trajectory_parallel

Protein contact networks
-------------------------

These functions underlie the ``protein-clustering`` command.

.. autofunction:: martinisoup.protein_clustering.build_nodes

.. autofunction:: martinisoup.protein_clustering.analyse_frame

.. autofunction:: martinisoup.protein_clustering.analyse_trajectory

.. autofunction:: martinisoup.protein_clustering.analyse_trajectory_parallel

Protein RDF
-----------

These functions underlie the ``protein-rdf`` command.

.. autofunction:: martinisoup.compute_rdf

.. autofunction:: martinisoup.compute_rdf_parallel
