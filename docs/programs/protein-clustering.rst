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
     --r-max FLOAT              Contact distance cutoff in Å (default: 6.0)
     --start INT                Start frame (default: 0)
     --stop INT                 Stop frame, exclusive (default: end of trajectory)
     --step INT                 Frame stride (default: 1)
     --output PATH              Output pickle file (default: protein_contacts.pkl)
     --parallel                 Run analysis in parallel
     --n_workers INT             Number of worker processes for parallel mode (default: 4)
     --chunk_size INT            Frames per chunk in parallel mode (default: 100)

Output pickle keys:

- ``command`` — the full command used to produce the results
- ``nodes`` — ``[{'id': int, 'type': str}, ...]``, one entry per protein molecule
- ``r_max`` — the contact distance cutoff used
- ``frames`` — list of per-frame NetworkX node-link graph dicts, each with ``frame``, ``time``,
  ``nodes``, and ``links`` keys
