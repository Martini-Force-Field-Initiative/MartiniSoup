"""
Metabolite–protein contact analysis using freud AABBQuery, built on
MDAnalysis's AnalysisBase (serial, multiprocessing or dask execution).

The input protein_types dictionary must have the form:

    {
        protein_name: {
            "segid": str,
            "n_monomers": int,
            "mask": np.ndarray(bool)
        },
        ...
    }

The mask must select atoms from the `proteins` AtomGroup belonging to
that protein type. The CLI prepares this dictionary.
"""

from __future__ import annotations
import numpy as np
import freud
from collections import defaultdict
from typing import Dict, Any, List

from MDAnalysis.analysis.base import ResultsGroup
from pint import Quantity

from .analysis_base import ParallelAnalysisBase, n_parts_for_chunk_size
from .units import unit_string


# --------------------------------------------------------------------------- #
# Utility: merge multiple partial results
# --------------------------------------------------------------------------- #

def _merge_counts(dict_list: List[Dict[str, Dict[str, int]]]) -> Dict[str, Dict[str, int]]:
    merged = defaultdict(lambda: defaultdict(int))
    for d in dict_list:
        for pname, resdict in d.items():
            for resname, count in resdict.items():
                merged[pname][resname] += count
    return {p: dict(r) for p, r in merged.items()}


# --------------------------------------------------------------------------- #
# Core per-frame analysis
# --------------------------------------------------------------------------- #

def _count_contacts_frame(
    proteins,
    metabolites,
    protein_types: Dict[str, Dict[str, Any]],
    cutoff_radius: float,
    metab_resnames: np.ndarray,
    metab_resids: np.ndarray,
    box_length: float,
) -> Dict[str, Dict[str, int]]:
    """Compute metabolite-protein contacts per protein type for one frame."""
    box = freud.Box.cube(box_length)
    aq = freud.locality.AABBQuery(box, metabolites.positions)

    frame_counts: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))

    for pname, info in protein_types.items():
        mask = info["mask"]

        prot_positions = proteins.positions[mask]
        if prot_positions.size == 0:
            continue

        result = aq.query(prot_positions, {"r_max": cutoff_radius})
        nlist = result.toNeighborList()

        if len(nlist.point_indices) == 0:
            continue

        met_indices = np.asarray(nlist.point_indices, dtype=np.int32)

        rids = metab_resids[met_indices]
        rnames = metab_resnames[met_indices]

        # Unique residues (resname, resid) pairs
        pairs = np.char.add(rnames.astype(str), "|" + rids.astype(str))
        uniq_pairs = np.unique(pairs)

        uniq_resnames = np.array([p.split("|", 1)[0] for p in uniq_pairs])

        names, freq = np.unique(uniq_resnames, return_counts=True)
        for resname, count in zip(names, freq):
            frame_counts[pname][resname] += int(count)

    return {p: dict(d) for p, d in frame_counts.items()}


# --------------------------------------------------------------------------- #
# AnalysisBase implementation (serial, multiprocessing, dask)
# --------------------------------------------------------------------------- #

class ContactCounts(ParallelAnalysisBase):
    """
    Metabolite-protein contact counts per protein type.

    For every frame, counts the unique metabolite residues within
    ``cutoff_radius`` of each protein type, and sums those counts over the
    trajectory. In a parallel run each part sums its own frames, the parts'
    counts are added together with :func:`_merge_counts`, and every protein
    type is guaranteed a (possibly empty) entry in the result.

    Parameters
    ----------
    proteins : MDAnalysis.AtomGroup
    metabolites : MDAnalysis.AtomGroup
    protein_types : dict
        Dictionary structured as::

            {
                protein_name: {
                    "segid": str,
                    "n_monomers": int,
                    "mask": boolean array (len == len(proteins))
                }
            }
    cutoff_radius : float or pint.Quantity
        Contact cutoff. A plain number is in the Universe's length unit (Å for
        a Universe read with default settings).
    **kwargs
        Passed to AnalysisBase (e.g. verbose).

    Results
    -------
    results.counts : dict
        {protein_name: {resname: count}}, one entry per protein type.
    results.cutoff_radius : float
        The cutoff used, in the Universe's length unit.
    results.units : dict
        {result name: pint unit string}; counts are dimensionless.

    Examples
    --------
    >>> cc = ContactCounts(proteins, metabolites, protein_types, cutoff_radius=5.0)
    >>> cc.run(backend='multiprocessing', n_workers=4)
    >>> cc.results.counts
    """

    def __init__(self, proteins, metabolites, protein_types, cutoff_radius, **kwargs):
        super().__init__(proteins.universe.trajectory, **kwargs)
        self._proteins = proteins
        self._metabolites = metabolites
        self._protein_types = protein_types
        self._cutoff_radius = self.unit_system.length_magnitude(cutoff_radius)
        self._metab_resnames = np.asarray(metabolites.atoms.resnames)
        self._metab_resids = np.asarray(metabolites.atoms.resids)

    def _prepare(self):
        # Plain dicts only: results are pickled back from worker processes.
        self.results.counts = {p: {} for p in self._protein_types}

    def _single_frame(self):
        frame_counts = _count_contacts_frame(
            self._proteins, self._metabolites, self._protein_types,
            self._cutoff_radius, self._metab_resnames, self._metab_resids,
            box_length=self._ts.dimensions[0],
        )
        counts = self.results.counts
        for pname, resdict in frame_counts.items():
            for resname, count in resdict.items():
                counts[pname][resname] = counts[pname].get(resname, 0) + count

    def _get_aggregator(self):
        return ResultsGroup(lookup={'counts': _merge_counts})

    def _conclude(self):
        # _merge_counts drops protein types with no contacts; restore them so
        # serial and parallel runs return the same keys.
        merged = self.results.counts
        self.results.counts = {p: merged.get(p, {}) for p in self._protein_types}
        self.results.cutoff_radius = self._cutoff_radius
        self.results.units = {
            'counts': 'dimensionless',
            'cutoff_radius': unit_string(self.unit_system.length_unit),
        }


# --------------------------------------------------------------------------- #
# Function wrappers (used by the binding-frequency command)
# --------------------------------------------------------------------------- #

def count_contacts(
    u,
    proteins,
    metabolites,
    protein_types: Dict[str, Dict[str, Any]],
    cutoff_radius: float | Quantity,
    start: int | None = None,
    stop: int | None = None,
    step: int = 1,
) -> Dict[str, Dict[str, int]]:
    """
    Compute metabolite-protein contacts per protein type (serial).

    Thin wrapper around :class:`ContactCounts` run with the serial backend.

    Parameters
    ----------
    u : MDAnalysis.Universe
        Unused beyond consistency with other analyses; the trajectory is taken
        from ``proteins.universe``.
    proteins : AtomGroup
    metabolites : AtomGroup
    protein_types : dict
        See :class:`ContactCounts`.
    cutoff_radius : float or pint.Quantity
        Contact cutoff. A plain number is in the Universe's length unit (Å
        for a Universe read with default settings).
    start, stop, step : frame iteration parameters

    Returns
    -------
    dict
        {protein_name: {resname: count}}
    """
    cc = ContactCounts(proteins, metabolites, protein_types, cutoff_radius, verbose=True)
    cc.run(start=start, stop=stop, step=step)
    return cc.results.counts


def count_contacts_parallel(
    topology: str,
    trajectory: str,
    prot_sel: str,
    metab_sel: str,
    protein_types: Dict[str, Dict[str, Any]],
    cutoff_radius: float | Quantity,
    start: int = 0,
    stop: int | None = None,
    step: int = 1,
    n_workers: int = 4,
    chunk_size: int = 100,
) -> Dict[str, Dict[str, int]]:
    """
    Compute metabolite-protein contacts per protein type in parallel.

    Thin wrapper around :class:`ContactCounts` run with MDAnalysis's
    ``multiprocessing`` backend. The frame range is split into contiguous
    parts of roughly ``chunk_size`` frames, distributed over ``n_workers``
    processes.

    Parameters
    ----------
    topology, trajectory : str
    prot_sel, metab_sel : selection strings
        ``protein_types`` masks must index into the ``prot_sel`` selection.
    protein_types : dict (see :class:`ContactCounts`)
    cutoff_radius : float or pint.Quantity
        Contact cutoff. A plain number is in the Universe's length unit (Å
        for a Universe read with default settings).
    start, stop, step : frame iteration parameters
    n_workers : int
    chunk_size : int
        Approximate frames per part.

    Returns
    -------
    dict
        {protein_name: {resname: count}}
    """
    import MDAnalysis as mda

    u = mda.Universe(topology, trajectory)
    n_parts = n_parts_for_chunk_size(u.trajectory, start, stop, step, chunk_size)

    cc = ContactCounts(
        u.select_atoms(prot_sel), u.select_atoms(metab_sel),
        protein_types, cutoff_radius,
    )
    cc.run(
        start=start, stop=stop, step=step,
        backend='multiprocessing', n_workers=n_workers, n_parts=n_parts,
    )
    return cc.results.counts
