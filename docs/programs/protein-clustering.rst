protein-clustering
===================

Builds a per-frame protein–protein contact network from a trajectory. Contacts are computed
at atom level and lifted to molecule level, so each unique protein-protein pair appears once
per frame. Output is a list of NetworkX node-link graphs, one per frame.

.. code-block:: bash

   martinisoup protein-clustering topology.tpr trajectory.xtc --output protein_contacts.pkl

.. code-block:: text

   positional arguments:
     topology              Topology file (e.g. .tpr, .gro)
     trajectory            Trajectory file (e.g. .xtc, .trr)

   options:
     --protein-selection STR   MDAnalysis selection string for proteins
     --r-max QUANTITY           Contact distance cutoff (default: 6 Å). A bare number is
                                in Å; a value with units, e.g. "0.6 nm", is converted.
     --start INT                Start frame (default: 0)
     --stop INT                 Stop frame, exclusive (default: end of trajectory)
     --step INT                 Frame stride (default: 1)
     --output PATH              Output pickle file (default: protein_contacts.pkl)
     --parallel                 Run analysis in parallel
     --n_workers INT             Number of worker processes for parallel mode (default: 4)
     --chunk_size INT            Approximate frames per part in parallel mode (default: 100)

Output pickle keys:

- ``command`` — the full command used to produce the results
- ``nodes`` — ``[{'id': int, 'type': str}, ...]``, one entry per protein molecule
- ``r_max`` — the contact distance cutoff used (Å)
- ``frames`` — list of per-frame NetworkX node-link graph dicts, each with ``frame``, ``time``
  (ps), ``nodes``, and ``links`` keys
- ``units`` — ``{key: pint unit string}``, with dotted paths for per-frame fields:
  ``r_max``, ``frames.frame``, ``frames.time``

Every output pickle includes a ``units`` dict giving the unit of each saved value; see :doc:`../units`.
