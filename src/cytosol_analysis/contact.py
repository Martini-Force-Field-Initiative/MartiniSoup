"""
Metabolite–protein contact analysis using freud AABBQuery,
with optional parallel execution and residue-level output.

The input protein_types dictionary must have the form:

    {
        protein_name: {
            "segid": str,
            "n_monomers": int,
            "mask": np.ndarray(bool)
        },
        ...
    }
"""

from __future__ import annotations
import numpy as np
import freud
from tqdm import tqdm
from collections import defaultdict
from typing import Dict, Any, List, Tuple


# --------------------------------------------------------------------------- #
# Utility: merge multiple partial results
# --------------------------------------------------------------------------- #

def _merge_atom_counts(dict_list):
    merged = defaultdict(lambda: defaultdict(int))
    for d in dict_list:
        for pname, resdict in d.items():
            for resname, count in resdict.items():
                merged[pname][resname] += count
    return {p: dict(r) for p, r in merged.items()}


def _merge_residue_counts(dict_list):
    merged = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    for d in dict_list:
        for pname, metdict in d.items():
            for met, resdict in metdict.items():
                for resid_key, count in resdict.items():
                    merged[pname][met][resid_key] += count
    return {
        p: {m: dict(r) for m, r in md.items()}
        for p, md in merged.items()
    }


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
    track_residues: bool = False,
) -> Dict[str, Any]:
    """
    Compute metabolite–protein contacts per protein type.

    Returns atom-level counts, and optionally residue-level contacts.
    """

    atom_results = {p: defaultdict(int) for p in protein_types.keys()}

    if track_residues:
        residue_results = {
            p: defaultdict(lambda: defaultdict(int))
            for p in protein_types.keys()
        }

    # Metabolite metadata (indexed by atom)
    metab_resnames = np.asarray(metabolites.atoms.resnames)
    metab_resids   = np.asarray(metabolites.atoms.resids)

    for ts in tqdm(u.trajectory[start:stop:step]):

        box = freud.Box.cube(u.dimensions[0])
        aq = freud.locality.AABBQuery(box, metabolites.positions)

        for pname, pdata in protein_types.items():

            mask = pdata["mask"]
            prot_atoms = proteins[mask]

            if len(prot_atoms) == 0:
                continue

            result = aq.query(prot_atoms.positions, {"r_max": cutoff_radius})
            nlist = result.toNeighborList()

            if len(nlist.point_indices) == 0:
                continue

            met_idx = np.asarray(nlist.point_indices, dtype=np.int32)

            # ------------------------------
            # Atom-level metabolite counting
            # ------------------------------
            met_resnames_frame = metab_resnames[met_idx]
            unique_met_resnames, counts = np.unique(
                met_resnames_frame, return_counts=True
            )

            for resname, count in zip(unique_met_resnames, counts):
                atom_results[pname][resname] += int(count)

            # ------------------------------
            # Residue-level protein contacts
            # ------------------------------
            if track_residues:
                prot_atoms_frame = prot_atoms[nlist.query_point_indices]

                # Build unique contacts per frame
                frame_contacts = set(
                    zip(
                        metab_resnames[met_idx],
                        prot_atoms_frame.resnames,
                        prot_atoms_frame.resids,
                    )
                )

                for met_res, prot_res, prot_id in frame_contacts:
                    key = (prot_res, int(prot_id))
                    residue_results[pname][met_res][key] += 1

    if track_residues:
        return {
            "atom_level": {p: dict(d) for p, d in atom_results.items()},
            "residue_level": {
                p: {m: dict(r) for m, r in md.items()}
                for p, md in residue_results.items()
            },
        }

    return {p: dict(d) for p, d in atom_results.items()}


# --------------------------------------------------------------------------- #
# Parallel version (multiprocessing)
# --------------------------------------------------------------------------- #

def _worker_contact_chunk(args):
    """
    Worker function for processing a chunk of frames.
    Rebuilds the Universe inside the worker.
    """
    (
        topology,
        trajectory,
        frame_range,
        prot_sel,
        metab_sel,
        protein_types,
        cutoff_radius,
        track_residues,
    ) = args

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
        track_residues=track_residues,
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
    track_residues: bool = False,
) -> Dict[str, Any]:
    """
    Parallel contact analysis by splitting trajectory into frame chunks.
    """

    import MDAnalysis as mda
    from concurrent.futures import ProcessPoolExecutor

    u = mda.Universe(topology, trajectory)
    n_frames = len(u.trajectory)

    frame_ids = np.arange(n_frames)
    chunks = [frame_ids[i:i + chunk_size]
              for i in range(0, n_frames, chunk_size)]

    worker_args = [
        (
            topology,
            trajectory,
            list(chunk),
            prot_sel,
            metab_sel,
            protein_types,
            cutoff_radius,
            track_residues,
        )
        for chunk in chunks
    ]

    atom_partials = []
    residue_partials = []

    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        for res in ex.map(_worker_contact_chunk, worker_args):
            if track_residues:
                atom_partials.append(res["atom_level"])
                residue_partials.append(res["residue_level"])
            else:
                atom_partials.append(res)

    atom_merged = _merge_atom_counts(atom_partials)

    if track_residues:
        residue_merged = _merge_residue_counts(residue_partials)
        return {
            "atom_level": atom_merged,
            "residue_level": residue_merged,
        }

    return atom_merged
