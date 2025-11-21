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
```python
from binding_analysis import TrajectoryAnalyzer, ResidenceAnalysis, KineticModels

# set up the groups
u = MDAnalysis.Universe(topology, trajectory)
metabolites = u.select_atoms("resname ATP")
proteins = u.select_atoms("protein")

analyzer = TrajectoryAnalyzer(u, metabolites, proteins, cutoff=4.0)
residences = analyzer.analyze()

analysis = ResidenceAnalysis(residences)
sc = analysis.survival_curve(moltype_id=0)
fit = KineticModels.fit_exponential(sc["time"], sc["survival"])
print(fit.fit_report())
```
