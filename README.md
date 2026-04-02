# MartiniSoup

Package for analysing metabolite interactions in heterogeneous cytosolic MD simulations.

## Features

- **Binding frequency analysis** — counts metabolite–protein contacts across a trajectory, normalised by protein monomer count and metabolite copy number
- **Residence time analysis** — tracks per-molecule binding events and extracts duration distributions for each metabolite type
- **Survival curve fitting** — empirical survival functions with bootstrap confidence intervals and single-exponential kinetic model fitting
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

Both programs are available as subcommands of the `martinisoup` CLI.

```
martinisoup <command> [options]

Commands:
  binding-frequency    Metabolite–protein contact frequency analysis
  residence-times      Binding event residence time analysis
```

### `binding-frequency`

Counts how often each metabolite type contacts each protein type over a trajectory. Results are normalised per protein monomer per metabolite copy.

```bash
martinisoup binding-frequency topology.tpr trajectory.xtc --pickle_out results.pkl
```

```
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
  --pickle_out PATH     Path to save results as a pickle file
```

The output pickle contains a dictionary with keys:
- `protein_results` — per-protein, per-metabolite raw and normalised contact counts
- `n_metabolites` — metabolite copy numbers in the system
- `n_proteins` — protein monomer counts in the system

### `residence-times`

Tracks binding and unbinding events for each metabolite molecule across the trajectory and records event durations.

```bash
martinisoup residence-times topology.tpr trajectory.xtc --output residences.pkl
```

```
positional arguments:
  topology              Topology file (e.g. .pdb, .gro, .tpr)
  trajectory            Trajectory file (e.g. .xtc, .dcd, .trr)

options:
  --metabolite-selection STR   MDAnalysis selection string for metabolites
  --protein-selection STR      MDAnalysis selection string for proteins
  --cutoff FLOAT               Distance cutoff in Å (default: 5.0)
  --start INT                  Start frame (default: 0)
  --stop INT                   Stop frame (default: end of trajectory)
  --step INT                   Frame step size (default: 1)
  --output PATH                Output pickle file (default: residues.pkl)
```

The output pickle contains a `{molecule_type: [durations]}` dictionary.

## Python API

```python
from martinisoup import BindingEventTracker, ResidenceAnalysis

# --- Residence time analysis ---
tracker = BindingEventTracker(u, metabolites, proteins, cutoff=5.0)
durations = tracker.track()  # {moltype_name: [durations]}

# --- Survival curve and fitting ---
ra = ResidenceAnalysis(durations)
t, S = ra.compute_survival("ATP")
lower, upper = ra.compute_survival_ci("ATP", n_boot=500)
result = ra.fit_exponential("ATP")

# --- Histogram ---
edges, centers, counts, pdf = ra.histogram("ATP", nbins=30, log_bins=True)
```
