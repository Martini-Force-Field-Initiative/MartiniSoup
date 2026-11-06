Python API
==========

The public API is importable directly from the ``martinisoup`` package:

.. code-block:: python

   from martinisoup import (
       ResidenceTracker, track_serial, track_parallel,
       compute_msd, BindingState, ResidenceRegistry,
       ProteinRDF, compute_rdf, compute_rdf_parallel,
   )

Parallel execution
-------------------

Every trajectory analysis in the package is an
``MDAnalysis.analysis.base.AnalysisBase`` subclass, built on
``martinisoup.analysis_base.ParallelAnalysisBase``:

- :class:`martinisoup.ProteinRDF` (``protein-rdf``)
- :class:`martinisoup.contact.ContactCounts` (``binding-frequency``)
- :class:`martinisoup.protein_clustering.ProteinContactNetwork` (``protein-clustering``)
- :class:`martinisoup.metabolite_partitioning.MetabolitePartitioning` (``metabolite-partitioning``)
- :class:`martinisoup.ResidenceTracker` (``residence-times``)

Each can be run directly on AtomGroups with any of the ``serial``,
``multiprocessing`` or ``dask`` backends:

.. code-block:: python

   rdf = ProteinRDF(u.select_atoms("protein"), r_max=100.0, n_bins=200)
   rdf.run(start=0, stop=None, step=1, backend="multiprocessing", n_workers=4)

The frame range is split into ``n_parts`` contiguous parts (default: one per
worker). Each part runs ``_prepare`` and ``_single_frame`` on its own copy of
the analysis, the parts' ``results`` are merged by the class's aggregator,
and ``_conclude`` runs once in the main process. The function-style
``*_parallel`` entry points used by the command-line programs are thin
wrappers. They build the Universe from file paths and convert ``chunk_size``
into ``n_parts`` (about ``chunk_size`` frames per part).

``ParallelAnalysisBase`` replaces MDAnalysis's stock ``multiprocessing``
backend with one that sends each part to a worker as its own task. The stock
backend (checked against MDAnalysis 2.10 and 2.11) batches several parts onto one
analysis object when ``n_parts > 4 * n_workers``, which silently returns
duplicated results. New analyses should subclass ``ParallelAnalysisBase``
rather than ``AnalysisBase`` directly.

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

``ResidenceTracker`` only parallelises the per-frame contact query (finding
metabolite atoms within cutoff of a protein atom). The binding/unbinding
state machine that turns those per-frame contacts into durations is
inherently sequential, because a binding event can span many frames and
therefore many parallel parts. It runs once, in the main process, in
``_conclude``, after the per-frame contact sets from all parts have been
concatenated in frame order, so serial and parallel runs give identical
durations. ``track_serial`` and ``track_parallel`` are thin wrappers that
return ``results.durations``.

.. autoclass:: martinisoup.ResidenceTracker

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

These underlie the ``binding-frequency`` command and can be called directly
for custom workflows. ``ContactCounts`` is an
``MDAnalysis.analysis.base.AnalysisBase`` subclass supporting the ``serial``,
``multiprocessing`` and ``dask`` backends; ``count_contacts`` and
``count_contacts_parallel`` are thin wrappers around it that return
``results.counts``. Every protein type in ``protein_types`` appears in the
result, with an empty dict if it made no contacts. ``count_contacts_parallel``
accepts the same
``start``/``stop``/``step`` frame-range parameters as ``count_contacts``
(previously it always processed the whole trajectory); the CLI does not yet
expose these for ``binding-frequency``.

.. autoclass:: martinisoup.contact.ContactCounts

.. autofunction:: martinisoup.contact.count_contacts

.. autofunction:: martinisoup.contact.count_contacts_parallel

.. autofunction:: martinisoup.scripts.binding_frequency.index_protein_segments

Database
--------

.. autofunction:: martinisoup.database.load_metabolite_classes

Metabolite partitioning
------------------------

These underlie the ``metabolite-partitioning`` command.
``MetabolitePartitioning`` is an ``MDAnalysis.analysis.base.AnalysisBase``
subclass supporting the ``serial``, ``multiprocessing`` and ``dask``
backends; ``analyse_trajectory`` and ``analyse_trajectory_parallel`` are thin
wrappers around it that return ``results.frames``.

.. autoclass:: martinisoup.metabolite_partitioning.MetabolitePartitioning

.. autofunction:: martinisoup.metabolite_partitioning.analyse_frame

.. autofunction:: martinisoup.metabolite_partitioning.analyse_trajectory

.. autofunction:: martinisoup.metabolite_partitioning.analyse_trajectory_parallel

Protein contact networks
-------------------------

These underlie the ``protein-clustering`` command.
``ProteinContactNetwork`` is an ``MDAnalysis.analysis.base.AnalysisBase``
subclass supporting the ``serial``, ``multiprocessing`` and ``dask``
backends; ``analyse_trajectory`` and ``analyse_trajectory_parallel`` are thin
wrappers around it that return ``results.frames``.

.. autoclass:: martinisoup.protein_clustering.ProteinContactNetwork

.. autofunction:: martinisoup.protein_clustering.build_nodes

.. autofunction:: martinisoup.protein_clustering.analyse_frame

.. autofunction:: martinisoup.protein_clustering.analyse_trajectory

.. autofunction:: martinisoup.protein_clustering.analyse_trajectory_parallel

Protein RDF
-----------

These underlie the ``protein-rdf`` command. ``ProteinRDF`` is an
``MDAnalysis.analysis.base.AnalysisBase`` subclass, so it supports the
``serial``, ``multiprocessing`` and ``dask`` backends directly:

.. code-block:: python

   import MDAnalysis as mda
   from martinisoup import ProteinRDF

   u = mda.Universe("topol.tpr", "traj.xtc")
   rdf = ProteinRDF(u.select_atoms("protein"), r_max=100.0, n_bins=200)
   rdf.run(backend="multiprocessing", n_workers=4)
   rdf.results.r, rdf.results.gr

``compute_rdf`` and ``compute_rdf_parallel`` are thin wrappers that return
the results as a plain dict.

.. autoclass:: martinisoup.ProteinRDF

.. autofunction:: martinisoup.compute_rdf

.. autofunction:: martinisoup.compute_rdf_parallel
