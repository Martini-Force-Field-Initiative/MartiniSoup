Python API
==========

The public API is importable directly from the ``martinisoup`` package:

.. code-block:: python

   from martinisoup import (
       track_serial, track_parallel, compute_msd, BindingState, ResidenceRegistry,
       compute_rdf, compute_rdf_parallel,
   )

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
directly for custom workflows.

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
