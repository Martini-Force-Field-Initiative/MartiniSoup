from collections import defaultdict

import freud
import numpy as np
from MDAnalysis.analysis.base import ResultsGroup
from MDAnalysis.lib.distances import capped_distance

from .analysis_base import ParallelAnalysisBase, n_parts_for_chunk_size


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


class MetabolitePartitioning(ParallelAnalysisBase):
    """
    Per-frame partitioning of metabolites into protein_adsorbed, clustered
    and soluble states.

    For every frame, each metabolite residue is classified by
    :func:`analyse_frame` and the per-resname counts are normalised to
    fractions. In a parallel run each part classifies its own frames and the
    parts' lists are concatenated in frame order.

    Parameters
    ----------
    metabolites : MDAnalysis.AtomGroup
        Metabolite atoms. Requires ``molnums`` (e.g. from a GROMACS .tpr).
    proteins : MDAnalysis.AtomGroup
    r_max : float
        Cluster neighbour cutoff in Angstrom (default: 5.0).
    contact_cutoff : float
        Protein contact cutoff in Angstrom (default: 5.0).
    **kwargs
        Passed to AnalysisBase (e.g. verbose).

    Results
    -------
    results.frames : list of dict
        One entry per frame, in frame order, with 'frame', 'time' and
        'fractions' keys; 'fractions' is
        {resname: {'protein_adsorbed': float, 'clustered': float, 'soluble': float}}.

    Examples
    --------
    >>> mp = MetabolitePartitioning(metabolites, proteins, r_max=5.0)
    >>> mp.run(backend='multiprocessing', n_workers=4)
    >>> mp.results.frames
    """

    def __init__(self, metabolites, proteins, r_max=5.0, contact_cutoff=5.0, **kwargs):
        super().__init__(metabolites.universe.trajectory, **kwargs)
        self._metabolites = metabolites
        self._proteins = proteins
        self._r_max = r_max
        self._contact_cutoff = contact_cutoff

    def _prepare(self):
        self.results.frames = []

    def _single_frame(self):
        counts = analyse_frame(
            self._metabolites,
            self._proteins,
            box_length=self._ts.dimensions[0],
            r_max=self._r_max,
            contact_cutoff=self._contact_cutoff,
        )
        # _normalize returns plain dicts; analyse_frame's defaultdict (built
        # with a lambda) could not be pickled back from a worker.
        self.results.frames.append({
            'frame': self._ts.frame,
            'time': self._ts.time,
            'fractions': _normalize(counts),
        })

    def _get_aggregator(self):
        return ResultsGroup(lookup={'frames': ResultsGroup.flatten_sequence})


def analyse_trajectory(u, metabolites, proteins,
                       r_max=5.0, contact_cutoff=5.0,
                       start=0, stop=None, step=1):
    """
    Run per-frame cluster state analysis over a trajectory slice (serial).

    Thin wrapper around :class:`MetabolitePartitioning` run with the serial
    backend.

    Parameters
    ----------
    u : MDAnalysis.Universe
        Unused beyond consistency with other analyses; the trajectory is taken
        from ``metabolites.universe``.
    metabolites : MDAnalysis.AtomGroup
    proteins : MDAnalysis.AtomGroup
    r_max : float
        Cluster neighbour cutoff in Angstrom.
    contact_cutoff : float
        Protein contact cutoff in Angstrom.
    start, stop, step : int
        Frame range, as in ``u.trajectory[start:stop:step]``.

    Returns
    -------
    list of dict
        Each entry has 'frame', 'time', and 'fractions' keys, in frame order.
    """
    mp = MetabolitePartitioning(
        metabolites, proteins, r_max=r_max, contact_cutoff=contact_cutoff,
        verbose=True,
    )
    mp.run(start=start, stop=stop, step=step)
    return mp.results.frames


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

    Thin wrapper around :class:`MetabolitePartitioning` run with MDAnalysis's
    ``multiprocessing`` backend. The frame range is split into contiguous
    parts of roughly ``chunk_size`` frames, distributed over ``n_workers``
    processes.

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
    start, stop, step : int
        Frame range, as in ``u.trajectory[start:stop:step]``.
    n_workers : int
        Number of worker processes.
    chunk_size : int
        Approximate frames per part.

    Returns
    -------
    list of dict
        Each entry has 'frame', 'time', and 'fractions' keys, in frame order.
    """
    import MDAnalysis as mda

    u = mda.Universe(topology, trajectory)
    n_parts = n_parts_for_chunk_size(u.trajectory, start, stop, step, chunk_size)

    mp = MetabolitePartitioning(
        u.select_atoms(metab_sel), u.select_atoms(prot_sel),
        r_max=r_max, contact_cutoff=contact_cutoff,
    )
    mp.run(
        start=start, stop=stop, step=step,
        backend='multiprocessing', n_workers=n_workers, n_parts=n_parts,
    )
    return mp.results.frames
