import pytest
import numpy as np
import MDAnalysis as mda
from martinisoup.contact import _merge_counts, count_contacts, count_contacts_parallel


def _make_universe(protein_positions, metabolite_positions,
                   met_resnames=None, met_resids=None, box_length=100.0):
    """
    Create a minimal single-frame MDAnalysis universe.

    Protein atoms get resname 'ALA'. Metabolite resnames and resids are
    configurable so that multi-atom residues can be tested.

    Parameters
    ----------
    protein_positions : list of (x, y, z)
    metabolite_positions : list of (x, y, z)
    met_resnames : list of str, one per metabolite atom (default: 'ATP' for all)
    met_resids : list of int, one per metabolite atom (default: sequential)
    box_length : float
        Side length of the cubic simulation box in Å.
    """
    n_prot = len(protein_positions)
    n_met = len(metabolite_positions)

    if met_resnames is None:
        met_resnames = ['ATP'] * n_met
    if met_resids is None:
        met_resids = list(range(1, n_met + 1))

    # Map unique (resname, resid) pairs → residue index, preserving insertion order
    met_residue_map = {}
    met_atom_resindex = []
    for rname, rid in zip(met_resnames, met_resids):
        key = (rname, rid)
        if key not in met_residue_map:
            met_residue_map[key] = len(met_residue_map)
        met_atom_resindex.append(n_prot + met_residue_map[key])

    n_met_residues = len(met_residue_map)
    n_residues = n_prot + n_met_residues

    u = mda.Universe.empty(
        n_prot + n_met,
        n_residues=n_residues,
        n_segments=2,
        atom_resindex=list(range(n_prot)) + met_atom_resindex,
        residue_segindex=[0] * n_prot + [1] * n_met_residues,
        trajectory=True,
    )

    # Per-residue attributes
    u.add_TopologyAttr('resname', ['ALA'] * n_prot + [k[0] for k in met_residue_map])
    u.add_TopologyAttr('resid', list(range(1, n_prot + 1)) + [k[1] for k in met_residue_map])

    # Per-atom attributes
    u.add_TopologyAttr('name', ['CA'] * n_prot + ['C1'] * n_met)

    u.atoms.positions = np.array(
        list(protein_positions) + list(metabolite_positions), dtype=float
    )
    u.dimensions = np.array([box_length, box_length, box_length, 90.0, 90.0, 90.0])
    return u


def _all_protein_types(proteins):
    """Return a protein_types dict treating all protein atoms as one type."""
    return {
        'PROT': {
            'mask': np.ones(len(proteins), dtype=bool),
            'n_monomers': 1,
        }
    }


class TestMergeCounts:

    def test_basic_additive_merge(self):
        a = {'PROT_A': {'ATP': 5, 'GTP': 2}}
        b = {'PROT_A': {'ATP': 3}}
        result = _merge_counts([a, b])
        assert result['PROT_A']['ATP'] == 8
        assert result['PROT_A']['GTP'] == 2

    def test_empty_list_returns_empty_dict(self):
        assert _merge_counts([]) == {}

    def test_single_dict_is_returned_unchanged(self):
        a = {'PROT_A': {'ATP': 4}}
        result = _merge_counts([a])
        assert result == {'PROT_A': {'ATP': 4}}

    def test_disjoint_proteins_are_combined(self):
        a = {'PROT_A': {'ATP': 3}}
        b = {'PROT_B': {'ATP': 7}}
        result = _merge_counts([a, b])
        assert result['PROT_A']['ATP'] == 3
        assert result['PROT_B']['ATP'] == 7

    def test_three_dicts_sum_correctly(self):
        dicts = [{'P': {'X': i}} for i in range(1, 4)]
        result = _merge_counts(dicts)
        assert result['P']['X'] == 6

    def test_disjoint_metabolites_in_same_protein(self):
        a = {'PROT_A': {'ATP': 2}}
        b = {'PROT_A': {'GTP': 5}}
        result = _merge_counts([a, b])
        assert result['PROT_A']['ATP'] == 2
        assert result['PROT_A']['GTP'] == 5


class TestCountContacts:

    def test_atom_within_cutoff_is_counted(self):
        u = _make_universe([[10.0, 10.0, 10.0]], [[12.0, 10.0, 10.0]])  # 2 Å apart
        proteins = u.select_atoms('resname ALA')
        metabolites = u.select_atoms('resname ATP')
        result = count_contacts(u, proteins, metabolites,
                                _all_protein_types(proteins), cutoff_radius=5.0)
        assert result['PROT']['ATP'] == 1

    def test_atom_outside_cutoff_is_not_counted(self):
        u = _make_universe([[10.0, 10.0, 10.0]], [[30.0, 30.0, 30.0]])  # ~35 Å apart
        proteins = u.select_atoms('resname ALA')
        metabolites = u.select_atoms('resname ATP')
        result = count_contacts(u, proteins, metabolites,
                                _all_protein_types(proteins), cutoff_radius=5.0)
        assert result['PROT'].get('ATP', 0) == 0

    def test_multiple_atoms_same_residue_count_once(self):
        """Two atoms from the same metabolite residue both within cutoff → 1 contact."""
        u = _make_universe(
            [[10.0, 10.0, 10.0]],
            [[11.0, 10.0, 10.0], [12.0, 10.0, 10.0]],
            met_resnames=['ATP', 'ATP'],
            met_resids=[1, 1],
        )
        proteins = u.select_atoms('resname ALA')
        metabolites = u.select_atoms('resname ATP')
        result = count_contacts(u, proteins, metabolites,
                                _all_protein_types(proteins), cutoff_radius=5.0)
        assert result['PROT']['ATP'] == 1

    def test_different_residues_count_separately(self):
        """Two atoms from different residues both within cutoff → 2 contacts."""
        u = _make_universe(
            [[10.0, 10.0, 10.0]],
            [[11.0, 10.0, 10.0], [13.0, 10.0, 10.0]],
            met_resnames=['ATP', 'ATP'],
            met_resids=[1, 2],
        )
        proteins = u.select_atoms('resname ALA')
        metabolites = u.select_atoms('resname ATP')
        result = count_contacts(u, proteins, metabolites,
                                _all_protein_types(proteins), cutoff_radius=5.0)
        assert result['PROT']['ATP'] == 2

    def test_empty_mask_gives_no_contacts(self):
        u = _make_universe([[10.0, 10.0, 10.0]], [[12.0, 10.0, 10.0]])
        proteins = u.select_atoms('resname ALA')
        metabolites = u.select_atoms('resname ATP')
        protein_types = {'PROT': {'mask': np.zeros(len(proteins), dtype=bool), 'n_monomers': 0}}
        result = count_contacts(u, proteins, metabolites,
                                protein_types, cutoff_radius=5.0)
        assert result['PROT'].get('ATP', 0) == 0

    def test_two_protein_types_counted_independently(self):
        """Two protein atoms with separate masks each get their own count."""
        u = _make_universe(
            [[10.0, 10.0, 10.0], [10.0, 15.0, 10.0]],
            [[12.0, 10.0, 10.0]],
        )
        proteins = u.select_atoms('resname ALA')
        metabolites = u.select_atoms('resname ATP')
        protein_types = {
            'TYPE_A': {'mask': np.array([True, False]), 'n_monomers': 1},
            'TYPE_B': {'mask': np.array([False, True]), 'n_monomers': 1},
        }
        result = count_contacts(u, proteins, metabolites,
                                protein_types, cutoff_radius=5.0)
        # Metabolite at (12,10,10): 2 Å from TYPE_A protein, 5 Å from TYPE_B protein
        assert result['TYPE_A']['ATP'] == 1
        assert result['TYPE_B'].get('ATP', 0) == 0


def _make_multiframe_universe(n_frames, n_prot=3, n_met=3, box_length=100.0, seed=0):
    """Multi-frame universe with random positions, for parallel/serial comparisons."""
    from MDAnalysis.coordinates.memory import MemoryReader

    rng = np.random.default_rng(seed)
    n_atoms = n_prot + n_met

    u = mda.Universe.empty(
        n_atoms, n_residues=n_atoms, n_segments=2,
        atom_resindex=list(range(n_atoms)),
        residue_segindex=[0] * n_prot + [1] * n_met,
        trajectory=True,
    )
    u.add_TopologyAttr('name', ['CA'] * n_prot + ['C1'] * n_met)
    u.add_TopologyAttr('resname', ['ALA'] * n_prot + ['ATP'] * n_met)
    u.add_TopologyAttr('resid', list(range(1, n_atoms + 1)))

    coords = rng.uniform(0, box_length, size=(n_frames, n_atoms, 3))
    dims = np.tile([box_length] * 3 + [90.0, 90.0, 90.0], (n_frames, 1))
    u.trajectory = MemoryReader(coords, order='fac', dimensions=dims, dt=1.0)
    return u


def _write_temp_files(u, tmp_path):
    gro = str(tmp_path / 'top.gro')
    xtc = str(tmp_path / 'traj.xtc')
    u.trajectory[0]
    u.atoms.write(gro)
    with mda.Writer(xtc, n_atoms=u.atoms.n_atoms) as w:
        for ts in u.trajectory:
            w.write(u.atoms)
    return gro, xtc


class TestCountContactsParallel:

    def test_matches_serial_across_uneven_chunk_boundaries(self, tmp_path):
        u = _make_multiframe_universe(37, n_prot=3, n_met=4)
        gro, xtc = _write_temp_files(u, tmp_path)

        u_direct = mda.Universe(gro, xtc)
        proteins = u_direct.select_atoms('resname ALA')
        metabolites = u_direct.select_atoms('resname ATP')
        protein_types = _all_protein_types(proteins)

        serial = count_contacts(u_direct, proteins, metabolites,
                                protein_types, cutoff_radius=20.0)
        parallel = count_contacts_parallel(
            gro, xtc, 'resname ALA', 'resname ATP', protein_types,
            cutoff_radius=20.0, n_workers=2, chunk_size=9,
        )
        assert serial == parallel

    def test_respects_start_stop_step(self, tmp_path):
        u = _make_multiframe_universe(20, n_prot=2, n_met=2)
        gro, xtc = _write_temp_files(u, tmp_path)

        u_direct = mda.Universe(gro, xtc)
        proteins = u_direct.select_atoms('resname ALA')
        metabolites = u_direct.select_atoms('resname ATP')
        protein_types = _all_protein_types(proteins)

        serial = count_contacts(u_direct, proteins, metabolites,
                                protein_types, cutoff_radius=20.0,
                                start=2, stop=18, step=3)
        parallel = count_contacts_parallel(
            gro, xtc, 'resname ALA', 'resname ATP', protein_types,
            cutoff_radius=20.0, start=2, stop=18, step=3,
            n_workers=2, chunk_size=3,
        )
        assert serial == parallel
