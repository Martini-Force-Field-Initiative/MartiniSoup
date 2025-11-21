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
    Compact tracker for one atom's binding events.

    Attributes
    ----------
    bound : bool
        True if currently bound.
    start : int | float | None
        Frame index or time when the current binding started.
    durations : list
        Completed durations (frame counts or time intervals).
    moltype_id : int | None
        Small integer id of the molecule type for this atom.
    """
    __slots__ = ("bound", "start", "durations", "moltype_id")
    def __init__(self, moltype_id=None):
        self.bound = False
        self.start = None
        self.durations = []
        self.moltype_id = moltype_id

    def update_bound(self, is_bound, current_stamp):
        """
        Update binding state at the current frame/time stamp.

        Parameters
        ----------
        is_bound : bool
            Whether the atom is considered bound at the current stamp.
        current_stamp : int or float
            Current frame index or time.
        """
        if is_bound:
            # start a binding event if not already bound
            if not self.bound:
                self.bound = True
                self.start = current_stamp
        else:
            # end a binding event if currently bound
            if self.bound:
                # store duration (end - start)
                self.durations.append(current_stamp - self.start)
                self.bound = False
                self.start = None

    def finalize(self, final_stamp):
        """
        Close any open binding event at the final stamp.

        Parameters
        ----------
        final_stamp : int or float
            Frame/time stamp to use as the unbinding time for still-bound atoms.
        """
        if self.bound:
            self.durations.append(final_stamp - self.start)
            self.bound = False
            self.start = None


class MetaboliteResidences:
    """
    Store residence-time data aggregated by atom, molecule, and molecule type.

    Parameters
    ----------
    tracker : dict[int, BoundState]
        Atom-level bound state tracking information.
    molnums : list[int]
        Molecule index for each atom (atom->molecule mapping).
    moltype_table : dict[int, str] or None
        Optional mapping from moltype_id to human-readable names.
    """

    def __init__(self, tracker, molnums, moltype_table=None, moltype_ids=None):
        self.tracker = tracker
        self.molnums = molnums  # per-atom → molecule mapping
        self.moltype_table = moltype_table  # int → type name
        self.moltype_ids = moltype_ids  # per-atom → type id

        self.molecule_data = self._aggregate_by_molecule()
        self.type_data = self._aggregate_by_type()

    # -------------------------------------------------------------
    # Aggregation helpers
    # -------------------------------------------------------------
    def _aggregate_by_molecule(self):
        """Aggregate durations by molecule index."""
        mol_dict = {}
        for atom_idx, state in self.tracker.items():
            mol = self.molnums[atom_idx]
            if mol not in mol_dict:
                mol_dict[mol] = []
            mol_dict[mol].extend(state.durations)
        return mol_dict

    def _aggregate_by_type(self):
        """Aggregate durations by moltype_id."""
        type_dict = {}
        for atom_idx, state in self.tracker.items():
            mt = state.moltype_id
            if mt not in type_dict:
                type_dict[mt] = []
            type_dict[mt].extend(state.durations)
        return type_dict

    def make_type_agg_named(self):
        """
        Create a dictionary mapping molecule-type names to aggregated durations.
        Requires `moltype_table` to be provided.
        """
        if self.moltype_table is None:
            raise ValueError("moltype_table is required for named aggregation")

        type_agg_named = {}
        for moltype_id, durations in self.type_data.items():
            name = self.moltype_table.get(moltype_id, f"type_{moltype_id}")
            if name not in type_agg_named:
                type_agg_named[name] = []
            type_agg_named[name].extend(durations)
        return type_agg_named

    # -------------------------------------------------------------
    # Serialization helpers for dictionary-only saving (pickle-safe)
    # -------------------------------------------------------------
    def to_dict(self):
        """Return a pure Python dictionary suitable for safe pickling."""
        out = {
            "molecule_data": self.molecule_data,
            "type_data": self.type_data,
            "molnums": self.molnums,
            "atoms": {},
            "moltype_table": self.moltype_table,
        }

        for idx, state in self.tracker.items():
            out["atoms"][idx] = {
                "durations": list(state.durations),
                "moltype_id": state.moltype_id,
            }

        return out

    @classmethod
    def from_dict(cls, d):
        """
        Reconstruct a lightweight MetaboliteResidences-like object
        from dictionary data. Does not reconstruct BoundState objects.
        """
        obj = cls.__new__(cls)
        obj.molecule_data = d["molecule_data"]
        obj.type_data = d["type_data"]
        obj.molnums = d["molnums"]
        obj.moltype_table = d.get("moltype_table", None)
        obj.tracker = None
        return obj
