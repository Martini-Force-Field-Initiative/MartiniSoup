"""
Fixtures for end-to-end tests that run the real ``martinisoup`` commands.

The test system is a small "cytosol" written to disk as a GROMACS ``.itp``
topology (which, unlike .gro, gives MDAnalysis molecule numbers and types, as
a .tpr would) plus ``.gro`` and ``.xtc`` coordinates. Residue names match the
CLIs' default protein and metabolite selections, so commands run without
selection options.
"""

import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use('Agg')

import numpy as np
import pytest
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader

from martinisoup.cli import main as cli_main

N_PROT, N_ATP, N_GLC = 4, 40, 40
N_FRAMES, DT_PS, BOX_A = 150, 10.0, 60.0
# Per-axis random-walk step of the metabolites, in Å per frame. Free 3D
# diffusion with this step has D = STEP_A**2 / (2 * DT_PS) Å²/ps.
STEP_A = 0.5
D_EXPECTED_CM2_S = STEP_A ** 2 / (2 * DT_PS) * 1e-4

_ITP = """\
[ defaults ]
1 1

[ atomtypes ]
P5  72.0  0.000  A  0.0  0.0
SN  54.0  0.000  A  0.0  0.0

[ moleculetype ]
PROT  1
[ atoms ]
1  P5  1  ALA  BB   1  0.0  72.0
2  P5  1  ALA  SC1  2  0.0  72.0
3  P5  2  GLY  BB   3  0.0  72.0

[ moleculetype ]
ATP  1
[ atoms ]
1  SN  1  ATP  C1  1  0.0  54.0
2  SN  1  ATP  C2  2  0.0  54.0

[ moleculetype ]
GLC  1
[ atoms ]
1  SN  1  GLC  C1  1  0.0  54.0

[ system ]
integration test cytosol

[ molecules ]
PROT  {n_prot}
ATP   {n_atp}
GLC   {n_glc}
"""


@dataclass(frozen=True)
class System:
    topology: str
    trajectory: str
    gro: str
    workdir: Path


def _coordinates(rng):
    """(n_frames, n_atoms, 3) positions in Å, wrapped into the box."""
    # Proteins: compact 3-bead bodies drifting slowly.
    prot_centres = rng.uniform(10, BOX_A - 10, (N_PROT, 3))
    prot_offsets = np.array([[0.0, 0, 0], [2.5, 0, 0], [0, 2.5, 0]])
    drift = np.cumsum(rng.normal(0, 0.1, (N_FRAMES, N_PROT, 1, 3)), axis=0)
    prot = (prot_centres[None, :, None, :] + prot_offsets[None, None] + drift)
    prot = prot.reshape(N_FRAMES, N_PROT * 3, 3)

    # Metabolites: free random walks with per-axis step STEP_A. Some start
    # close to proteins so that contacts and residence events occur.
    n_met = N_ATP + N_GLC
    near = prot_centres[rng.integers(0, N_PROT, n_met // 2)] + rng.normal(0, 3, (n_met // 2, 3))
    far = rng.uniform(0, BOX_A, (n_met - n_met // 2, 3))
    starts = np.vstack([near, far])
    walk = starts + np.cumsum(rng.normal(0, STEP_A, (N_FRAMES, n_met, 3)), axis=0)
    atp = np.stack([walk[:, :N_ATP], walk[:, :N_ATP] + [2.0, 0, 0]], axis=2)
    atp = atp.reshape(N_FRAMES, N_ATP * 2, 3)
    glc = walk[:, N_ATP:]

    return np.concatenate([prot, atp, glc], axis=1) % BOX_A


@pytest.fixture(scope='session')
def system(tmp_path_factory):
    d = tmp_path_factory.mktemp('cytosol')
    itp = d / 'system.itp'
    itp.write_text(_ITP.format(n_prot=N_PROT, n_atp=N_ATP, n_glc=N_GLC))

    u = mda.Universe(str(itp))
    coords = _coordinates(np.random.default_rng(42))
    assert coords.shape[1] == u.atoms.n_atoms
    dims = np.tile([BOX_A] * 3 + [90.0] * 3, (N_FRAMES, 1))
    u.load_new(coords, format=MemoryReader, order='fac', dimensions=dims, dt=DT_PS)

    gro, xtc = d / 'conf.gro', d / 'traj.xtc'
    u.atoms.write(str(gro))
    with mda.Writer(str(xtc), n_atoms=u.atoms.n_atoms) as w:
        for ts in u.trajectory:
            w.write(u.atoms)
    return System(topology=str(itp), trajectory=str(xtc), gro=str(gro), workdir=d)


@pytest.fixture
def run_cli(monkeypatch):
    """Run ``martinisoup <command> <args...>`` in-process, as from the shell."""
    def run(*args):
        monkeypatch.setattr(sys, 'argv', ['martinisoup', *map(str, args)])
        cli_main()
    return run
