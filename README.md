# MartiniSoup

Package for analysing metabolite interactions in heterogeneous cytosolic MD simulations.

## Features

- **Binding frequency analysis** — counts metabolite–protein contacts across a trajectory, normalised by protein monomer count and metabolite copy number
- **Residence time analysis** — tracks per-molecule binding events and extracts duration distributions for each metabolite type
- **Mean squared displacement** — computes MSD per metabolite type using the Einstein relation, averaged across all molecules of each type
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
  msd                  Mean squared displacement per metabolite type
  msd-fitter           Fit MSD data to extract diffusion coefficients
```

By default, the MDAnalysis selection strings for metabolites and proteins are:

```
        metabolite_selection = 'not resname NA CL ION GLU ASN VAL ALA GLY ARG SER PRO THR PHE GLN LYS LEU ASP ILE MET HIS CYS TRP TYR'

        protein_selection = 'resname GLU ASN VAL ALA GLY ARG SER PRO THR PHE GLN LYS LEU ASP ILE MET HIS CYS TRP TYR'

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

### `msd`

Computes the mean squared displacement per metabolite type using the MDAnalysis Einstein MSD implementation. A NoJump transformation is applied automatically to unwrap periodic boundary conditions. Results are averaged across all copies of each molecule type.

```bash
martinisoup msd topology.tpr trajectory.xtc --output msd.pkl
```

```
positional arguments:
  topology              Topology file (e.g. .tpr, .gro)
  trajectory            Trajectory file (e.g. .xtc, .trr)

options:
  --metab_sel STR       MDAnalysis selection string for metabolites
  --start INT           Start frame (default: first frame of second half of trajectory)
  --stop INT            Stop frame (default: start + 100)
  --fft                 Use FFT algorithm for MSD calculation (default: True)
  --output PATH         Output pickle file (default: msd.pkl)
```

The output pickle contains a dictionary with keys:
- `resnames` — list of unique metabolite type names
- `residue_timeseries` — mean MSD per type, shape `(n_types, n_lagtimes)`
- `residue_std` — standard deviation across molecules, shape `(n_types, n_lagtimes)`
- `time` — absolute times corresponding to each lagtime
- `lagtimes` — MSD values averaged across all molecules
- `dimensions` — dimensionality factor used in the MSD calculation
- `dt` — trajectory timestep in ps

### `msd-fitter`

Fits MSD output from `martinisoup msd` to a linear model and extracts diffusion coefficients via the Einstein relation (D = slope / 6). Multiple input files are treated as replicas and averaged before fitting.

```bash
# single replica
martinisoup msd-fitter msd.pkl --output-dir results/

# average over replicas before fitting
martinisoup msd-fitter replica_1/msd.pkl replica_2/msd.pkl replica_3/msd.pkl --output-dir results/

# custom fitting window and matplotlib style
martinisoup msd-fitter msd.pkl --cut-start 10 --cut-end 50 --style mystyle.mplstyle --output-dir results/
```

```
positional arguments:
  files                 Pickle file(s) from `martinisoup msd`. Multiple files are averaged as replicas.

options:
  --cut-start INT       Start index of the fitting window (default: 10)
  --cut-end INT         End index of the fitting window (default: 50)
  --output-dir PATH     Directory for output plots and results CSV (default: current directory)
  --style PATH          Path to a matplotlib style file
```

Outputs per-residue MSD plots and a `diffusion_coefficients.csv` with columns `resname`, `D`, `D_err`.

When multiple replica files are provided, the mean MSD across replicas is computed and the uncertainty is propagated as `sqrt(sum(σ_i²)) / n_replicas` before fitting.

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
