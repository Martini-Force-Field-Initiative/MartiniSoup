protein-rdf
============

Computes the protein centre-of-mass radial distribution function over a trajectory.
Each protein molecule's centre of mass is computed per frame, then all pairwise distances
are histogrammed and normalised by the mean box volume and number of protein pairs.

.. code-block:: bash

   martinisoup protein-rdf topology.tpr trajectory.xtc --output protein_rdf.pkl

.. code-block:: text

   positional arguments:
     topology              Topology file (e.g. .tpr, .gro)
     trajectory            Trajectory file (e.g. .xtc, .trr)

   options:
     --protein-selection STR   MDAnalysis selection string for protein atoms
     --r-max FLOAT              Upper distance limit in Å (default: 100.0)
     --n-bins INT                Number of histogram bins (default: 200)
     --start INT                Start frame (default: 0)
     --stop INT                 Stop frame, exclusive (default: end of trajectory)
     --step INT                 Frame stride (default: 1)
     --output PATH              Output pickle file (default: protein_rdf.pkl)
     --parallel                 Run analysis in parallel
     --n_workers INT             Number of worker processes for parallel mode (default: 4)
     --chunk_size INT            Frames per chunk in parallel mode (default: 100)

Output pickle keys:

- ``command`` — the full command used to produce the results
- ``r`` — bin midpoints in Å
- ``gr`` — g(r) values
- ``n_proteins`` — number of protein molecules
- ``mean_vol`` — mean box volume over analysed frames in Å³
