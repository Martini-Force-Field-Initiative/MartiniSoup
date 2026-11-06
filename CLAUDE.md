# MartiniSoup

## Development

- Set up: `uv sync` (editable install plus dev tools). Run tests with
  `uv run pytest`.
- Tests marked `integration` (`tests/integration/`) run the real `martinisoup`
  commands on a small generated system written to disk: a GROMACS `.itp`
  topology (used instead of a `.tpr`, which MDAnalysis can't write; it still
  provides `molnums`/`moltypes`) with `.gro` and `.xtc` coordinates. Run only
  those with `uv run pytest -m integration`, or skip them with
  `-m "not integration"`.
- Docs: `uv run --with-requirements docs/requirements.txt sphinx-build -W
  --keep-going -b html docs docs/_build/html`. Warnings are errors, as in CI.
- CI (`.github/workflows/`): `tests.yml` runs unit and integration tests on
  Python 3.10–3.14 for pushes to `master` and pull requests. `docs.yml` builds
  the docs on pull requests and uploads the HTML as an artifact.
  `astral-sh/setup-uv` publishes no major tags from v8 on, so it is pinned to
  an exact release (e.g. `@v10.3.3`); `@v10` does not exist.
- Parallel workers are always started with `spawn`
  (`analysis_base.START_METHOD`), never `fork`. Forked workers inherit open
  MDAnalysis readers and threads; this segfaulted in CI (Linux, Python 3.10),
  and forcing `fork` on macOS deadlocks. `tests/test_analysis_base.py` guards
  this. Any script that starts parallel work needs an
  `if __name__ == "__main__":` guard, including throwaway debugging scripts.
- macOS + Python 3.14: if `uv run` cannot import `martinisoup`, the editable
  install's `.pth` file has the macOS "hidden" flag, which Python 3.14 skips.
  Clear it with `chflags nohidden .venv/lib/python3.14/site-packages/__editable__.martinisoup-*.pth`,
  or run with `PYTHONPATH=src`.

## Known bugs

### binding-frequency only works with .tpr segment names

**Status:** known, not yet fixed. `TestBindingFrequency` in
`tests/integration/test_cli.py` is a strict `xfail` for this bug. Once the bug
is fixed that test will pass, which a strict xfail reports as a failure, so
remove the marker then.

`index_protein_segments` (`src/martinisoup/scripts/binding_frequency.py`)
derives protein names from segment IDs, assuming the `.tpr` pattern
`seg_<n>_<name>`:

- Segment IDs without that pattern (e.g. from a GROMACS `.itp`/`.top`,
  where the segid is the molecule type, such as `PROT`) take the `else`
  branch, `protein_name = parts`. That is a *list*, and it is then used as a
  dict key, so the command crashes with `TypeError: unhashable type: 'list'`.
- Even with a `.tpr`, a protein type listed in more than one block of
  `[ molecules ]` gets one segment per block (`seg_0_X`, `seg_3_X`). Both map
  to name `X`, and the later block overwrites the earlier one's mask and
  monomer count, so contacts and normalisation silently cover only the last
  block.

Likely fix: group protein atoms by `moltypes` (available from both `.tpr` and
`.itp`), combining all molecules of a type into one mask, with `n_monomers` =
number of distinct `molnums` of that type. This changes how names are derived
for `.tpr` input (moltype instead of segid suffix; normally identical), so
note it in `docs/changelog.rst`.

### residence-fitter crashes on histograms with empty bins

`fit_and_plot` (`src/martinisoup/scripts/residence_fitter.py`) fits a
`PowerLawModel` to the binned counts. If any bin in the `--bins-*` range is
empty, the fit fails (`ValueError: The model function generated NaN values`),
because of zero counts on log axes. With log-spaced bins this easily happens
at long lifetimes or for rarely seen metabolites. Possible fixes: drop empty
bins before fitting (and from the weights), or fit in log space with zeros
excluded, and report which groups had too few populated bins to fit.

### Residence times do not track which protein a metabolite is bound to

**Status:** known, not yet fixed. Affects `residence-times`, `ResidenceTracker`,
`track_serial` and `track_parallel`, and every lifetime distribution or
`residence-fitter` exponent computed from them.

**Symptom.** A residence event is meant to be the time a metabolite stays bound
to *one* protein. What is measured is the time a metabolite stays continuously
in contact with *any* protein. All of these are recorded as a single unbroken
event:

1. Near protein A at frame *t*, near protein B at frame *t+1*. The switch is
   visible in the saved frames, but it is still merged.
2. In contact with A and B at the same time during a handover, then B only.
3. Leaves A and binds B between two saved frames.

**Cause.**

- `_bound_metabolite_atoms` (`src/martinisoup/residence_tracker.py`) builds one
  KDTree from *all* protein atoms and keeps only the metabolite side of each
  contact pair: `frozenset(pair[1] for pair in sdm.keys())`. The protein atom
  index (`pair[0]`), and so the protein's identity, is thrown away.
- `_run_state_machine` and `BindingState.update_bound`
  (`src/martinisoup/data_structures.py`) keep a single bound/unbound flag per
  metabolite atom, with no protein attached.

**Related problem: tracking is per atom (bead), not per molecule.** There is one
`BindingState` per metabolite *atom*, and
`ResidenceRegistry.get_durations_by_type` pools every atom's events by molecule
type. As a result:

- A multi-bead metabolite that binds once contributes one event per bead in
  contact.
- When one bead's contact flickers while the molecule stays bound, that bead's
  events split into many short ones, which biases the distribution towards
  short lifetimes.

Fix this together with the protein identity problem, because both change how
an event is defined.

**What can be fixed.** Cases 1 and 2 above, plus the per-bead counting:

- Map each contact's protein atom to its protein molecule (`molnums`) and each
  metabolite atom to its metabolite molecule. Track state per
  *(metabolite molecule, protein molecule)* pair.
- A pair counts as bound in a frame if any bead of the metabolite is within the
  cutoff of any atom of that protein. Its event ends at the first frame with no
  contact to that protein, even if the metabolite is touching another protein.
  A handover therefore gives two overlapping events, one per protein.
- Store state only for pairs currently in contact (a dict keyed by pair), not a
  dense metabolites × proteins array.
- The per-frame results passed through the parallel aggregation become sets of
  `(metabolite_mol, protein_mol)` pairs instead of sets of atom indices. The
  ordered fold in `_conclude` stays sequential, as now.

**What cannot be fixed from the trajectory.** Case 3 (and briefly unbinding from
a protein then rebinding to the same protein between frames) cannot be told
apart from continuous binding. Residence times are only resolved to the frame
interval `dt`; events or gaps shorter than `dt` are invisible. Document this
limitation. The only real remedy is writing frames more often. At most, an
event could be *flagged* as suspicious (for example, if the metabolite's
nearest protein changes during it), but what happened between frames can't be
recovered.

**Constraints and decisions for the fix.**

- Protein identity needs `molnums` (e.g. from a `.tpr`). The current fallback
  in `_molnums` (residue indices) would make every protein *residue* a separate
  "protein". When protein `molnums` are missing, raise an error; don't fall
  back silently.
- The fix changes the residence-time distributions. Keep the current "any
  protein, per atom" definition available behind an option (e.g.
  `--any-protein` on the CLI), so results can be compared with earlier
  analyses. Record which definition was used in the output (as `stamp` is
  recorded now).
- Output: consider also saving which protein each event was on, so per-protein
  analyses are possible. Keep the existing `{moltype: [durations]}` layout and
  the `units` / `stamp` / `dt` keys so that `residence-fitter` keeps working.
- Tests to add: a metabolite moving A → B in consecutive frames gives two
  events; an A+B handover gives two overlapping events; a multi-bead metabolite
  binding once gives one event; a single bead flickering inside a bound
  molecule does not split the event; serial and multiprocessing runs agree,
  including events that cross part boundaries (the existing
  `test_multiprocessing_matches_serial` pattern in
  `tests/test_residence_tracker.py`).
- Docs: update `docs/programs/residence.rst` (event definition, the `dt`
  resolution limit, the new option), and add an entry to `docs/changelog.rst`
  saying that distributions change compared with earlier versions.
