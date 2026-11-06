import numpy as np
import pytest
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.lib.distances import self_distance_array

from martinisoup.protein_rdf import ProteinRDF, compute_rdf


def _make_protein_universe(n_frames, n_mols=8, atoms_per_mol=3, box_length=80.0, seed=0):
    """In-memory universe of `n_mols` small 'proteins' with molnums and masses."""
    rng = np.random.default_rng(seed)
    n_atoms = n_mols * atoms_per_mol
    u = mda.Universe.empty(
        n_atoms, n_residues=n_mols, n_segments=1,
        atom_resindex=np.repeat(np.arange(n_mols), atoms_per_mol),
        residue_segindex=[0] * n_mols,
        trajectory=True,
    )
    u.add_TopologyAttr('name', ['BB'] * n_atoms)
    u.add_TopologyAttr('resname', ['ALA'] * n_mols)
    u.add_TopologyAttr('resid', list(range(1, n_mols + 1)))
    u.add_TopologyAttr('molnums', np.arange(n_mols))  # per-residue attr
    u.add_TopologyAttr('masses', rng.uniform(50, 100, size=n_atoms))

    # Molecule centres anywhere in the box (incl. across boundaries), atoms
    # jittered around them, box size fluctuating slightly between frames.
    centres = rng.uniform(0, box_length, size=(n_frames, n_mols, 1, 3))
    jitter = rng.normal(0, 1.5, size=(n_frames, n_mols, atoms_per_mol, 3))
    coords = (centres + jitter).reshape(n_frames, n_atoms, 3)
    lengths = box_length + rng.normal(0, 0.5, size=(n_frames, 1))
    dims = np.hstack([np.repeat(lengths, 3, axis=1), np.full((n_frames, 3), 90.0)])
    u.trajectory = MemoryReader(coords, order='fac', dimensions=dims, dt=1.0)
    return u


def _reference_rdf(u, atoms, r_max, n_bins, start, stop, step):
    """Independent brute-force implementation of the same normalisation."""
    mols = atoms.split('molecule')
    n = len(mols)
    bins = np.linspace(0, r_max, n_bins + 1)
    counts = np.zeros(n_bins)
    vols = []
    for ts in u.trajectory[start:stop:step]:
        coms = np.array([m.center_of_mass(wrap=True) for m in mols])
        counts += np.histogram(self_distance_array(coms, box=ts.dimensions), bins=bins)[0]
        vols.append(np.prod(ts.dimensions[:3]))
    mean_vol = np.mean(vols)
    shell = (4 / 3) * np.pi * (bins[1:] ** 3 - bins[:-1] ** 3)
    gr = (2 * mean_vol) / (n * (n - 1)) * counts / (shell * len(vols))
    return 0.5 * (bins[:-1] + bins[1:]), gr, mean_vol


class TestProteinRDF:

    def test_serial_matches_reference(self):
        u = _make_protein_universe(20)
        atoms = u.atoms
        rdf = ProteinRDF(atoms, r_max=40.0, n_bins=40).run(start=2, stop=18, step=3)

        r_ref, gr_ref, vol_ref = _reference_rdf(u, atoms, 40.0, 40, 2, 18, 3)
        assert np.allclose(rdf.results.r, r_ref)
        assert np.allclose(rdf.results.gr, gr_ref)
        assert rdf.results.mean_vol == pytest.approx(vol_ref)
        assert rdf.results.n_proteins == 8
        assert len(rdf.results.volumes) == len(range(2, 18, 3))

    @pytest.mark.parametrize("n_workers,n_parts", [(2, 2), (3, 7), (2, 23), (2, 24)])
    def test_multiprocessing_matches_serial(self, n_workers, n_parts):
        """
        Uneven part boundaries must not change counts, volumes or g(r).

        n_parts > 4 * n_workers is the regression case for the stock
        MDAnalysis multiprocessing backend, which batches several parts onto
        one analysis object and returns duplicated results (see
        martinisoup.analysis_base).
        """
        u = _make_protein_universe(53)
        serial = ProteinRDF(u.atoms, r_max=40.0, n_bins=40).run(start=3, stop=50, step=2)
        parallel = ProteinRDF(u.atoms, r_max=40.0, n_bins=40).run(
            start=3, stop=50, step=2,
            backend='multiprocessing', n_workers=n_workers, n_parts=n_parts,
        )

        assert np.array_equal(parallel.results.counts, serial.results.counts)
        assert np.allclose(parallel.results.volumes, serial.results.volumes)
        assert np.allclose(parallel.results.gr, serial.results.gr)
        assert parallel.results.mean_vol == pytest.approx(serial.results.mean_vol)

    def test_compute_rdf_wrapper_returns_expected_dict(self):
        u = _make_protein_universe(10)
        out = compute_rdf(u, u.atoms, r_max=40.0, n_bins=40)
        r_ref, gr_ref, vol_ref = _reference_rdf(u, u.atoms, 40.0, 40, 0, None, 1)

        assert set(out) == {'r', 'gr', 'n_proteins', 'mean_vol', 'units'}
        assert out['units'] == {
            'r': 'angstrom', 'gr': 'dimensionless',
            'n_proteins': 'dimensionless', 'mean_vol': 'angstrom ** 3',
        }
        assert np.allclose(out['r'], r_ref)
        assert np.allclose(out['gr'], gr_ref)
        assert out['mean_vol'] == pytest.approx(vol_ref)
        assert out['n_proteins'] == 8
