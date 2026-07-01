CLI Reference
=============

All commands are available as subcommands of the ``martinisoup`` CLI:

.. code-block:: text

   martinisoup <command> [options]

   Commands:
     binding-frequency        Metabolite–protein contact frequency analysis
     residence-times          Binding event residence time analysis
     msd                      Mean squared displacement per metabolite type
     msd-fitter               Fit MSD data to extract diffusion coefficients
     residence-fitter         Fit residence-time histograms by metabolite class
     metabolite-partitioning  Classify metabolites as clustered, soluble, or protein-adsorbed
     protein-clustering       Build per-frame protein–protein contact networks
     protein-rdf              Protein centre-of-mass radial distribution function

Default selections
------------------

Unless overridden, MDAnalysis selection strings default to:

.. code-block:: text

   metabolite_selection = 'not resname NA CL ION GLU ASN VAL ALA GLY ARG SER PRO THR PHE GLN LYS LEU ASP ILE MET HIS CYS TRP TYR'
   protein_selection    = 'resname GLU ASN VAL ALA GLY ARG SER PRO THR PHE GLN LYS LEU ASP ILE MET HIS CYS TRP TYR'

.. toctree::
   :maxdepth: 1
   :caption: Programs

   programs/binding-frequency
   programs/metabolite-partitioning
   programs/protein-clustering
   programs/protein-rdf
   programs/msd
   programs/residence
