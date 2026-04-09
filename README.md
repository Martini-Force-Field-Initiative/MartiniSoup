# MartiniSoup

Package for analysing metabolite interactions in heterogeneous cytosolic MD simulations.

Full documentation is available at [martinisoup.readthedocs.io](https://martinisoup.readthedocs.io).

## Features

- **Binding frequency analysis** — counts metabolite–protein contacts across a trajectory, normalised by protein monomer count and metabolite copy number
- **Residence time analysis** — tracks per-molecule binding events and extracts duration distributions for each metabolite type
- **Mean squared displacement** — computes MSD per metabolite type using the Einstein relation, averaged across all molecules of each type
- **Parallel trajectory processing** — optional multiprocessing support for large trajectories

## Installation

```bash
pip install .
```

For development:
```bash
pip install -e .
```

## Usage

Commands are available as subcommands of the `martinisoup` CLI. See the [CLI reference](https://martinisoup.readthedocs.io/en/latest/cli.html) for full option documentation.

```
martinisoup <command> [options]

Commands:
  binding-frequency    Metabolite–protein contact frequency analysis
  residence-times      Binding event residence time analysis
  msd                  Mean squared displacement per metabolite type
  msd-fitter           Fit MSD data to extract diffusion coefficients
  residence-fitter     Fit residence-time histograms by metabolite class
```

## Python API

```python
from martinisoup import track_serial, track_parallel, compute_msd

# --- Serial residence time tracking ---
residences = track_serial(
    topology, trajectory,
    metabolite_selection, protein_selection,
    cutoff=5.0, start=0, stop=None, step=1
)
# residences: {moltype_name: [durations]}

# --- Parallel residence time tracking ---
residences = track_parallel(
    topology, trajectory,
    metabolite_selection, protein_selection,
    cutoff=5.0, start=0, stop=None, step=1,
    n_workers=4, chunk_size=100
)

# --- MSD analysis ---
results = compute_msd(topology, trajectory, metabolite_selection)
```

See the [API reference](https://martinisoup.readthedocs.io/en/latest/api.html) for full documentation.
