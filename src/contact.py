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
# Core analysis (single process)
# --------------------------------------------------------------------------- #

def run_contact_analysis(
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
        Dictionary of:
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

        L = u.dimensions[0]
        box = freud.Box.cube(L)

        # Build spatial index for metabolites
        aq = freud.locality.AABBQuery(box, metabolites.positions)

        # Query once per protein type
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

            # Extract metadata
            rids = metab_resids[met_indices]
            rnames = metab_resnames[met_indices]

            # Unique residues (resname, resid) pairs
            pairs = np.char.add(rnames.astype(str), "|" + rids.astype(str))
            uniq_pairs = np.unique(pairs)

            uniq_resnames = np.array([p.split("|", 1)[0] for p in uniq_pairs])

            # Count per resname
            names, freq = np.unique(uniq_resnames, return_counts=True)

            for resname, count in zip(names, freq):
                results[pname][resname] += int(count)

    # Convert defaultdicts to dict
    return {p: dict(d) for p, d in results.items()}


# --------------------------------------------------------------------------- #
# Parallel version (multiprocessing)
# --------------------------------------------------------------------------- #

def _worker_contact_chunk(args):
    """
    Worker function for processing a chunk of frames.
    Rebuilds the Universe inside the worker.
    """
    (topology, trajectory, frame_range,
     prot_sel, metab_sel, protein_types, cutoff_radius) = args

    import MDAnalysis as mda

    u = mda.Universe(topology, trajectory)
    proteins = u.select_atoms(prot_sel)
    metabolites = u.select_atoms(metab_sel)

    return run_contact_analysis(
        u,
        proteins,
        metabolites,
        protein_types,
        cutoff_radius,
        start=frame_range[0],
        stop=frame_range[-1] + 1,
        step=1,
    )


def parallel_contact_analysis(
    topology: str,
    trajectory: str,
    prot_sel: str,
    metab_sel: str,
    protein_types: Dict[str, Dict[str, Any]],
    cutoff_radius: float,
    n_workers: int = 4,
    chunk_size: int = 100,
) -> Dict[str, Dict[str, int]]:
    """
    Parallel contact analysis by splitting trajectory into frame chunks.

    Parameters
    ----------
    topology, trajectory : str
    prot_sel, metab_sel : selection strings
    protein_types : dict (see above)
    cutoff_radius : float
    n_workers : int
    chunk_size : int

    Returns
    -------
    dict
        Merged contact counts.
    """

    import MDAnalysis as mda
    import numpy as np
    from concurrent.futures import ProcessPoolExecutor

    u = mda.Universe(topology, trajectory)
    n_frames = len(u.trajectory)

    frame_ids = np.arange(n_frames)
    chunks = [frame_ids[i:i+chunk_size] for i in range(0, n_frames, chunk_size)]

    worker_args = [
        (topology, trajectory, list(chunk),
         prot_sel, metab_sel, protein_types, cutoff_radius)
        for chunk in chunks
    ]

    partial_results = []

    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        for res in ex.map(_worker_contact_chunk, worker_args):
            partial_results.append(res)

    return _merge_counts(partial_results)
