"""
Metabolite–protein contact analysis using freud AABBQuery,
with optional parallel execution and metadata-rich output.

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
from tqdm import tqdm
from collections import defaultdict
from typing import Dict, Any, List

from .parallel import map_trajectory_parallel


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
# Core per-frame analysis (shared by serial and parallel paths)
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
# Core analysis (single process)
# --------------------------------------------------------------------------- #

def count_contacts(
    u,
    proteins,
    metabolites,
    protein_types: Dict[str, Dict[str, Any]],
    cutoff_radius: float,
    start: int | None = None,
    stop: int | None = None,
    step: int = 1,
) -> Dict[str, Dict[str, int]]:
    """
    Compute metabolite–protein contacts per protein type.

    Parameters
    ----------
    u : MDAnalysis.Universe
    proteins : AtomGroup
    metabolites : AtomGroup
    protein_types : dict
        Dictionary structured as::

            {
                protein_name: {
                    "segid": str,
                    "n_monomers": int,
                    "mask": boolean array (len == len(proteins))
                }
            }
    cutoff_radius : float
    start, stop, step : frame iteration parameters

    Returns
    -------
    dict
        {protein_name: {resname: count}}
    """

    results = {p: defaultdict(int) for p in protein_types.keys()}

    # Atom-level metabolite metadata
    metab_resnames = np.asarray(metabolites.atoms.resnames)
    metab_resids   = np.asarray(metabolites.atoms.resids)

    for ts in tqdm(u.trajectory[start:stop:step]):
        frame_counts = _count_contacts_frame(
            proteins, metabolites, protein_types, cutoff_radius,
            metab_resnames, metab_resids, box_length=u.dimensions[0],
        )
        for pname, resdict in frame_counts.items():
            for resname, count in resdict.items():
                results[pname][resname] += count

    # Convert defaultdicts to dict
    return {p: dict(d) for p, d in results.items()}


# --------------------------------------------------------------------------- #
# Parallel version (multiprocessing)
# --------------------------------------------------------------------------- #

def _contact_setup(u, prot_sel, metab_sel, protein_types, cutoff_radius):
    """Build the per-worker context: protein/metabolite atoms and metadata."""
    proteins = u.select_atoms(prot_sel)
    metabolites = u.select_atoms(metab_sel)
    metab_resnames = np.asarray(metabolites.atoms.resnames)
    metab_resids = np.asarray(metabolites.atoms.resids)
    return proteins, metabolites, protein_types, cutoff_radius, metab_resnames, metab_resids


def _contact_per_frame(u, context, ts):
    """Compute one frame's contact counts."""
    proteins, metabolites, protein_types, cutoff_radius, metab_resnames, metab_resids = context
    return _count_contacts_frame(
        proteins, metabolites, protein_types, cutoff_radius,
        metab_resnames, metab_resids, box_length=u.dimensions[0],
    )


def count_contacts_parallel(
    topology: str,
    trajectory: str,
    prot_sel: str,
    metab_sel: str,
    protein_types: Dict[str, Dict[str, Any]],
    cutoff_radius: float,
    start: int = 0,
    stop: int | None = None,
    step: int = 1,
    n_workers: int = 4,
    chunk_size: int = 100,
) -> Dict[str, Dict[str, int]]:
    """
    Parallel contact analysis by splitting the trajectory into frame chunks.

    Parameters
    ----------
    topology, trajectory : str
    prot_sel, metab_sel : selection strings
    protein_types : dict (see above)
    cutoff_radius : float
    start, stop, step : frame iteration parameters
    n_workers : int
    chunk_size : int

    Returns
    -------
    dict
        Merged contact counts.
    """
    per_frame_results = map_trajectory_parallel(
        topology, trajectory, _contact_setup, _contact_per_frame,
        setup_args=(prot_sel, metab_sel, protein_types, cutoff_radius),
        start=start, stop=stop, step=step,
        n_workers=n_workers, chunk_size=chunk_size,
        desc="Counting contacts (parallel)",
    )

    return _merge_counts(per_frame_results)
