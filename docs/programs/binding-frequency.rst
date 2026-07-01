binding-frequency
==================

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
