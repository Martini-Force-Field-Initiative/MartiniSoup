Changelog
=========

Unreleased
----------

These changes come from the units overhaul (see :doc:`units`). Several change
the contents of output files; the **Changed output** entries matter if you
have scripts that read MartiniSoup pickles directly.

Changed output
^^^^^^^^^^^^^^

- ``msd``: ``lagtimes`` now holds the lag time of each MSD point (ps). It
  previously held the MSD averaged across all particles (Å²), despite its
  name; that is now under the new ``msd`` key. ``msd-fitter`` still reads
  files in the old layout and converts them automatically. Scripts that read
  ``lagtimes`` from ``msd`` output should switch to ``msd``.
- ``msd``: the ``trajectory_units`` key has been removed. It described the
  trajectory file's units (e.g. nm for ``.xtc``), not the units of the saved
  data, which have always been Å and ps.
- ``msd-fitter``: ``units`` in ``diffusion_coefficients.pkl`` is now a dict,
  ``{'results': 'centimeter ** 2 / second'}``, instead of the string
  ``"cm²/s"``, matching every other output file.
- ``compute_rdf`` / ``compute_rdf_parallel`` (Python API) return an extra
  ``units`` key.

New
^^^

- Every output pickle includes a ``units`` dict giving the unit of each saved
  value as a pint-parseable string. Programs whose output had no unitful
  values now also save the cutoffs they used (``cutoff``, ``r_max``,
  ``contact_cutoff``).
- Distance, temperature and viscosity options accept values with units, e.g.
  ``--cutoff "0.5 nm"`` or ``--temperature "37 degC"``. A bare number keeps
  its previous meaning (Å, K or mPa·s), so existing commands are unchanged.
- Analysis classes accept pint Quantities for distances, record their output
  units in ``results.units``, and work correctly on Universes read with
  ``convert_units=False``.
- ``residence-times --frames`` measures durations in trajectory frames
  instead of time. ``residence-times`` output now also records ``stamp``
  (``'time'`` or ``'frame'``) and ``dt`` (time between frames).
- ``residence-fitter`` converts lifetimes from the units recorded in each
  file, including frame counts (using the saved ``dt``), instead of assuming
  ps.
- New module ``martinisoup.units`` with the shared pint registry and helpers;
  ``pint`` is a new dependency.

Fixed
^^^^^

- Parallel runs (``--parallel``, ``backend="multiprocessing"``) could crash
  with a segmentation fault, or hang, on Linux with Python 3.13 or earlier,
  where worker processes were started by forking. Workers are now always
  started with ``spawn``. Python scripts that run a parallel analysis must do
  so under an ``if __name__ == "__main__":`` guard on every platform (this was
  already required on macOS and Windows).

- ``msd-fitter`` interpreted MSDs using the trajectory file's units (e.g. nm
  for ``.xtc``) although the data was in Å, giving diffusion coefficients
  100 times too large for ``.xtc`` trajectories. In addition: fit weights
  are now 1/σ (previously 1/σ²); the Einstein factor is 2·d from the saved
  MSD dimensionality (previously always 6); and the Yeh-Hummer correction is
  only applied when ``--temperature`` is given (previously it defaulted to
  300 K).

Compatibility
^^^^^^^^^^^^^

Python 3.10 is no longer supported; MartiniSoup now requires Python 3.11 or
later.

Files written by earlier versions have no ``units`` entry. ``msd-fitter`` and
``residence-fitter`` still read them, with a warning, assuming Å and ps
(the units MartiniSoup has always written).
