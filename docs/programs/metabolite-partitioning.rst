metabolite-partitioning
=======================

Classifies each metabolite residue as ``protein_adsorbed``, ``clustered``, or ``soluble`` at
every frame of a trajectory. Clustering uses a neighbour-cutoff approach (via ``freud``) and
protein contacts are detected with a distance cutoff. A residue is classified as
``protein_adsorbed`` only if it has strictly more protein contacts than metabolite–metabolite
contacts; otherwise it falls back to the cluster-size rule (``soluble`` if its cluster contains
only that molecule, ``clustered`` otherwise).

.. code-block:: bash

   martinisoup metabolite-partitioning topology.tpr trajectory.xtc --output cluster_states.pkl

.. code-block:: text

   positional arguments:
     topology              Topology file (e.g. .tpr, .gro)
     trajectory            Trajectory file (e.g. .xtc, .trr)

   options:
     --metabolite-selection STR   MDAnalysis selection string for metabolites
     --protein-selection STR      MDAnalysis selection string for proteins
     --r-max FLOAT                Cluster neighbour cutoff in Å (default: 5.0)
     --contact-cutoff FLOAT       Protein contact cutoff in Å (default: 5.0)
     --start INT                  Start frame (default: 0)
     --stop INT                   Stop frame, exclusive (default: end of trajectory)
     --step INT                   Frame stride (default: 1)
     --output PATH                Output pickle file (default: cluster_states.pkl)
     --parallel                   Run analysis in parallel
     --n_workers INT              Number of worker processes for parallel mode (default: 4)
     --chunk_size INT             Frames per chunk in parallel mode (default: 100)

Output pickle keys:

- ``command`` — the full command used to produce the results
- ``results`` — list of per-frame dicts, each with ``frame``, ``time``, and ``fractions`` keys.
  ``fractions`` is ``{resname: {'protein_adsorbed': float, 'clustered': float, 'soluble': float}}``.
