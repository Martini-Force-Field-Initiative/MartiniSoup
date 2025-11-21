
"""
Data structures for storing and organizing residence time information.

This module defines:
    • BoundState – compact storage for atom-level binding events
    • MetaboliteResidences – aggregation of residence durations at
      atom, molecule, and molecule-type levels
"""

from collections import defaultdict

class BoundState:
    """
    Track binding state and residence times for a single atom.

    Attributes
    ----------
    bound : bool
        Whether the atom is currently bound.
    start : float or None
        Timestep when the current binding event started.
    durations : list of float
        Completed binding event durations.
    moltype_id : int or None
        ID describing the molecular type this atom belongs to.
    """

    __slots__ = ("bound", "start", "durations", "moltype_id")

    def __init__(self, moltype_id=None):
        self.bound = False
        self.start = None
        self.durations = []
        self.moltype_id = moltype_id


class MetaboliteResidences:
    """
    Container storing residence durations at the atom, molecule, and
    molecule-type level.

    Parameters
    ----------
    tracker : dict[int, BoundState]
        Dictionary keyed by atom index storing atom-level events.
    molnums : array-like
        Mapping from atom index to molecule index.
    moltype_table : dict[int, str]
        Mapping from moltype_id → molecule type name.

    Attributes
    ----------
    tracker : dict[int, BoundState]
    molecule_data : dict[int, dict]
        Aggregated residence info per molecule.
    type_data : dict[int, list[float]]
        Aggregated durations per molecule type.
    """

    def __init__(self, tracker, molnums, moltype_table):
        self.tracker = tracker
        self.molnums = molnums
        self.moltype_table = moltype_table

        self.molecule_data = self._aggregate_by_molecule()
        self.type_data = self._aggregate_by_moltype()

    # ------------------------------------------------------------------
    def _aggregate_by_molecule(self):
        """Aggregate atom-level durations into molecule-level durations."""

        mol_agg = defaultdict(lambda: {"moltype_id": None, "durations": []})

        for atom_idx, state in self.tracker.items():
            molnum = self.molnums[atom_idx]
            mt_id = state.moltype_id

            if mol_agg[molnum]["moltype_id"] is None:
                mol_agg[molnum]["moltype_id"] = mt_id

            mol_agg[molnum]["durations"].extend(state.durations)

        return dict(mol_agg)

    # ------------------------------------------------------------------
    def _aggregate_by_moltype(self):
        """Aggregate molecule-level durations into molecule-type-level lists."""

        type_agg = defaultdict(list)

        for molnum, entry in self.molecule_data.items():
            mt_id = entry["moltype_id"]
            type_agg[mt_id].extend(entry["durations"])

        return dict(type_agg)
