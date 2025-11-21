# binding_analysis

Tools for analyzing ligand–protein residence times from MD simulations.

## Features
- Efficient binding-event tracking across large trajectories
- Histograms (linear and log-binned)
- Survival curves and hazard estimation
- LMfit-accelerated kinetic model fitting
- Bootstrap confidence intervals

## Installation
```bash
pip install .
```

## Usage Example

### Command line example
A dataset containing all binding event times per metabolite can be generated with 
the following command:

```commandline
binding-analysis \
    topology.pdb trajectory.xtc \
    --metabolite-selection "resname ATP ADP AMP" \
    --protein-selection "protein" \
    --cutoff 5.0 \
    --output results.p
```

The results can then be investigated further by opening them with a script:

```python
import pickle

with open("results.p", "rb") as f:
    type_agg_named = pickle.load(f)

# Example:
print(list(type_agg_named.keys()))   # ['ATP', 'ADP', 'AMP']
for molecule, n_events in type_agg_named.items():
    print(f"Number of {molecule} events:", len(n_events))
# Number of ATP events: 10
# ...
```
