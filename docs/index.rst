MartiniSoup
===========

Tools for analysing metabolite interactions in heterogeneous cytosolic molecular dynamics simulations.

MartiniSoup provides a command-line interface and Python API for:

- **Binding frequency analysis** — counts metabolite–protein contacts across a trajectory, normalised by protein monomer count and metabolite copy number
- **Residence time analysis** — tracks per-molecule binding events and extracts duration distributions
- **Mean squared displacement** — computes MSD per metabolite type using the Einstein relation
- **Metabolite partitioning** — classifies metabolites as clustered, soluble, or protein-adsorbed at every frame
- **Protein contact networks** — builds per-frame protein–protein contact graphs for clustering analysis
- **Protein RDF** — computes the protein centre-of-mass radial distribution function
- **Parallel trajectory processing** — optional multiprocessing support for large trajectories
- **Data fitting** — power-law fitting for residence times and linear fitting for MSD

.. toctree::
   :maxdepth: 2
   :caption: Contents

   installation
   cli
   api
