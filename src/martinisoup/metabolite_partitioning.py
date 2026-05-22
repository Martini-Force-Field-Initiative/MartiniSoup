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

    # Union-Find: merge any clusters bridged by the same molecule.
    # Taking the minimum ID per molecule (the old approach) orphaned molecules
    # that were connected through a bridging multi-bead molecule.
    parent = {int(c): int(c) for c in np.unique(cluster_idx)}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]  # path compression
            x = parent[x]
        return x

    for mol in np.unique(molnums):
        mol_clusters = [int(c) for c in np.unique(cluster_idx[molnums == mol])]
        root = find(mol_clusters[0])
        for cid in mol_clusters[1:]:
            r = find(cid)
            if r != root:
                parent[max(root, r)] = min(root, r)
                root = min(root, r)

    cluster_idx = np.array([find(int(c)) for c in cluster_idx])

    # Cluster sizes in molecules (not atoms)
    cluster_to_mols = defaultdict(set)
    for cid, mol in zip(cluster_idx, molnums):
        cluster_to_mols[cid].add(mol)
    cluster_sizes = {cid: len(mols) for cid, mols in cluster_to_mols.items()}

    # Protein contacts per atom
    prot_pairs, _ = capped_distance(
        metabolites.positions,
        proteins.positions,
        max_cutoff=contact_cutoff,
        box=metabolites.dimensions,
    )
    prot_contacts = np.bincount(prot_pairs[:, 0], minlength=len(metabolites))

    # Metabolite-metabolite contacts per atom (cross-residue only)
    metab_pairs, _ = capped_distance(
        metabolites.positions,
        metabolites.positions,
        max_cutoff=r_max,
        box=metabolites.dimensions,
    )
    res_indices = metabolites.resindices
    cross = res_indices[metab_pairs[:, 0]] != res_indices[metab_pairs[:, 1]]
    metab_contacts = np.bincount(metab_pairs[cross, 0], minlength=len(metabolites))

    # Per-residue classification: protein_adsorbed only wins if it has strictly
    # more bead contacts than met-met contacts; ties fall through to cluster rule.
    counts = defaultdict(lambda: {'protein_adsorbed': 0, 'soluble': 0, 'clustered': 0})
    resnames = metabolites.resnames
    seen_residues = set()

    for i in range(len(metabolites)):
        ri = res_indices[i]
        if ri in seen_residues:
            continue
        seen_residues.add(ri)

        mask = res_indices == ri
        if prot_contacts[mask].sum() > metab_contacts[mask].sum():
            counts[resnames[i]]['protein_adsorbed'] += 1
        else:
            cid = cluster_idx[i]
            if cluster_sizes[cid] == 1:
                counts[resnames[i]]['soluble'] += 1
            else:
                counts[resnames[i]]['clustered'] += 1

    return counts


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


# --------------------------------------------------------------------------- #
# Parallel version (multiprocessing)
# --------------------------------------------------------------------------- #

def _worker_clustering_chunk(args):
    """
    Worker function for processing a chunk of frames.
    Rebuilds the Universe inside the worker process.
    """
    (topology, trajectory, frame_indices,
     metab_sel, prot_sel, r_max, contact_cutoff) = args

    import MDAnalysis as mda

    u = mda.Universe(topology, trajectory)
    metabolites = u.select_atoms(metab_sel)
    proteins = u.select_atoms(prot_sel)

    results = []
    for frame_idx in frame_indices:
        u.trajectory[frame_idx]
        counts = analyse_frame(
            metabolites,
            proteins,
            box_length=u.dimensions[0],
            r_max=r_max,
            contact_cutoff=contact_cutoff,
        )
        results.append({
            'frame': u.trajectory.ts.frame,
            'time': u.trajectory.ts.time,
            'fractions': _normalize(counts),
        })
    return results


def analyse_trajectory_parallel(
    topology,
    trajectory,
    metab_sel,
    prot_sel,
    r_max=5.0,
    contact_cutoff=5.0,
    start=0,
    stop=None,
    step=1,
    n_workers=4,
    chunk_size=100,
):
    """
    Parallel per-frame cluster state analysis.

    Each frame is classified independently, so the trajectory slice is split
    into chunks and dispatched to worker processes. Results are collected and
    returned sorted by frame index.

    Parameters
    ----------
    topology, trajectory : str
        Paths to topology and trajectory files.
    metab_sel, prot_sel : str
        MDAnalysis selection strings for metabolites and proteins.
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
    n_workers : int
        Number of worker processes.
    chunk_size : int
        Frames per chunk.

    Returns
    -------
    list of dict
        Each entry has 'frame', 'time', and 'fractions' keys, sorted by frame.
    """
    import MDAnalysis as mda
    from concurrent.futures import ProcessPoolExecutor

    u = mda.Universe(topology, trajectory)
    n_frames = len(u.trajectory)
    stop = stop if stop is not None else n_frames

    frame_indices = list(range(start, stop, step))
    chunks = [frame_indices[i:i + chunk_size]
              for i in range(0, len(frame_indices), chunk_size)]

    worker_args = [
        (topology, trajectory, chunk, metab_sel, prot_sel, r_max, contact_cutoff)
        for chunk in chunks
    ]

    all_results = []
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        for chunk_results in tqdm(executor.map(_worker_clustering_chunk, worker_args),
                                  total=len(chunks), desc="Analysing trajectory (parallel)"):
            all_results.extend(chunk_results)

    all_results.sort(key=lambda x: x['frame'])
    return all_results
