# MartiniSoup

Package for analysing metabolite interactions in heterogeneous cytosolic MD simulations.

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

Both programs are available as subcommands of the `martinisoup` CLI.

```
martinisoup <command> [options]

Commands:
  binding-frequency    Metabolite–protein contact frequency analysis
  residence-times      Binding event residence time analysis
  msd                  Mean squared displacement per metabolite type
```

By default, the MDAnalysis selection strings for metabolites and proteins are:

```
metabolite_selection = 'not resname NA CL ION GLU ASN VAL ALA GLY ARG SER PRO THR PHE GLN LYS LEU ASP ILE MET HIS CYS TRP TYR'

protein_selection = 'resname GLU ASN VAL ALA GLY ARG SER PRO THR PHE GLN LYS LEU ASP ILE MET HIS CYS TRP TYR'
```

### `binding-frequency`

Counts how often each metabolite type contacts each protein type over a trajectory. Results are normalised per protein monomer per metabolite copy.

```bash
martinisoup binding-frequency topology.tpr trajectory.xtc --output results.pkl
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
  --output PATH         Path to save results as a pickle file (default: binding.pkl)
```

The output pickle contains a dictionary with keys:
- `command` — the full command used to produce the results
- `protein_results` — per-protein, per-metabolite raw and normalised contact counts
- `n_metabolites` — metabolite copy numbers in the system
- `n_proteins` — protein monomer counts in the system

### `residence-times`

Tracks binding and unbinding events for each metabolite molecule across the trajectory and records event durations. The use of the `--summary` flag is strongly recommended, the raw data output may be quite large for data analysis purposes.

```bash
martinisoup residence-times topology.tpr trajectory.xtc --output lifetimes.pkl --summary summary.pkl
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
  --output PATH                Output pickle file (default: lifetimes.pkl)
  --summary FILE               Write a summary pickle (unique residence times and
                               their counts) to FILE
  --parallel                   Run contact detection in parallel (Phase 1), then
                               run the state machine sequentially (Phase 2).
                               Produces identical results to serial mode.
  --n_workers INT              Number of worker processes for parallel mode (default: 4)
  --chunk_size INT             Frames per chunk in parallel mode (default: 100)
```

The output pickle contains a dictionary with keys:
- `command` — the full command used to produce the results
- `residences` — `{molecule_type: [durations]}` dictionary

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
- `command` — the full command used to produce the results
- `resnames` — list of unique metabolite type names
- `residue_timeseries` — mean MSD per type, shape `(n_types, n_lagtimes)`
- `residue_std` — standard deviation across molecules, shape `(n_types, n_lagtimes)`
- `time` — absolute times corresponding to each lagtime
- `lagtimes` — MSD values averaged across all molecules
- `dimensions` — dimensionality factor used in the MSD calculation
- `dt` — trajectory timestep in ps

## Python API

```python
from martinisoup import track_serial, track_parallel

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
```
