from collections import defaultdict

import freud
import numpy as np
from MDAnalysis.lib.distances import capped_distance
from tqdm import tqdm


def analyse_frame(metabolites, proteins, box_length,
                  r_max=5.0, contact_cutoff=5.0):
    """
    Classify each metabolite residue as protein_adsorbed, clustered, or soluble
    for a single trajectory frame.

    Parameters
    ----------
    metabolites : MDAnalysis.AtomGroup
    proteins : MDAnalysis.AtomGroup
    box_length : float
        Cubic box side length in Angstrom.
    r_max : float
        Cluster neighbour cutoff in Angstrom.
    contact_cutoff : float
        Protein contact cutoff in Angstrom.

    Returns
    -------
    dict
        {resname: {'protein_adsorbed': int, 'clustered': int, 'soluble': int}}
    """
    # Clustering (atom-level)
    box = freud.Box.cube(box_length)
    cl = freud.cluster.Cluster()
    cl.compute((box, metabolites.positions), neighbors={'r_max': r_max})

    cluster_idx = cl.cluster_idx.copy()
    molnums = metabolites.molnums

    # Ensure one cluster ID per molecule (use minimum)
    for mol in np.unique(molnums):
        mask = molnums == mol
        cluster_idx[mask] = cluster_idx[mask].min()

    # Cluster sizes in molecules (not atoms)
    cluster_to_mols = defaultdict(set)
    for cid, mol in zip(cluster_idx, molnums):
        cluster_to_mols[cid].add(mol)
    cluster_sizes = {cid: len(mols) for cid, mols in cluster_to_mols.items()}

    # Protein contacts (vectorised)
    pairs, _ = capped_distance(
        metabolites.positions,
        proteins.positions,
        max_cutoff=contact_cutoff,
        box=metabolites.dimensions,
    )
    protein_bound_atoms = np.zeros(len(metabolites), dtype=bool)
    protein_bound_atoms[pairs[:, 0]] = True

    # Per-residue classification
    results = defaultdict(lambda: {'protein_adsorbed': 0, 'soluble': 0, 'clustered': 0})
    res_indices = metabolites.resindices
    resnames = metabolites.resnames
    seen_residues = set()

    for i in range(len(metabolites)):
        ri = res_indices[i]
        if ri in seen_residues:
            continue
        seen_residues.add(ri)

        resname = resnames[i]
        if protein_bound_atoms[res_indices == ri].any():
            results[resname]['protein_adsorbed'] += 1
        else:
            cid = cluster_idx[i]
            if cluster_sizes[cid] == 1:
                results[resname]['soluble'] += 1
            else:
                results[resname]['clustered'] += 1

    return results


def _normalize(counts):
    """Convert raw per-resname counts to fractions."""
    out = {}
    for k, v in counts.items():
        tot = sum(v.values())
        out[k] = {kk: vv / tot for kk, vv in v.items()}
    return out


def analyse_trajectory(u, metabolites, proteins,
                       r_max=5.0, contact_cutoff=5.0,
                       start=0, stop=None, step=1):
    """
    Run per-frame cluster state analysis over a trajectory slice.

    Parameters
    ----------
    u : MDAnalysis.Universe
    metabolites : MDAnalysis.AtomGroup
    proteins : MDAnalysis.AtomGroup
    r_max : float
        Cluster neighbour cutoff in Angstrom.
    contact_cutoff : float
        Protein contact cutoff in Angstrom.
    start : int
        First frame index.
    stop : int or None
        Last frame index (exclusive). None means end of trajectory.
    step : int
        Frame stride.

    Returns
    -------
    list of dict
        Each entry has 'frame', 'time', and 'fractions' keys.
    """
    results = []
    for ts in tqdm(u.trajectory[start:stop:step]):
        counts = analyse_frame(
            metabolites,
            proteins,
            box_length=u.dimensions[0],
            r_max=r_max,
            contact_cutoff=contact_cutoff,
        )
        results.append({
            'frame': ts.frame,
            'time': ts.time,
            'fractions': _normalize(counts),
        })
    return results
