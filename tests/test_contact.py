import pytest
import numpy as np
import MDAnalysis as mda
from martinisoup.contact import _merge_counts, count_contacts


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
