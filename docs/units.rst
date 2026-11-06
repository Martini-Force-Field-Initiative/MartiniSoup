Units
=====

MartiniSoup uses `pint <https://pint.readthedocs.io>`_ to track physical units.
Every program works out the units of the trajectory data it reads, does its
calculations in those units, and records the unit of every saved number, so
that later analysis, inside or outside MartiniSoup, can use the values
correctly.

Units of the trajectory data
----------------------------

MDAnalysis converts coordinates, box dimensions and times to its own units
(**Å** and **ps**) as it reads a trajectory. This happens whatever units the
file uses (``.xtc`` files store nm, for example). MartiniSoup therefore reads
the units of the data *as held in memory*:

- With a Universe read normally (the default, and what every
  ``martinisoup`` command does), lengths are in Å and times in ps.
- Only with a Universe created with ``convert_units=False`` is the data left
  in the file's own units (e.g. nm). MartiniSoup detects this and works in,
  and reports, those units.

.. note::

   ``u.trajectory.units`` describes the trajectory *file*, not the data in
   memory. It reports ``nm`` for an ``.xtc`` file even though the positions
   MDAnalysis returns are in Å. Using it to interpret the data was the
   cause of an earlier unit bug in ``msd-fitter``. Use
   :meth:`martinisoup.units.UnitSystem.from_universe` instead.

Giving values with units
------------------------

**On the command line**, options that take a distance, temperature or
viscosity accept either a bare number or a value with units:

.. code-block:: bash

   martinisoup binding-frequency topol.tpr traj.xtc --cutoff 5         # 5 Å
   martinisoup binding-frequency topol.tpr traj.xtc --cutoff "0.5 nm"  # also 5 Å
   martinisoup msd-fitter msd.pkl --temperature "37 degC"              # 310.15 K

A bare number is always in the unit stated in that option's help and in the
program's documentation:

=========================================================  =========================
Option                                                     Bare number means
=========================================================  =========================
``--cutoff``, ``--r-max``, ``--contact-cutoff``             Å
``--temperature`` (``msd-fitter``)                          K
``--eta`` (``msd-fitter``)                                  mPa·s
``--bins-start``, ``--bins-stop`` (``residence-fitter``)    log10 of the lifetime in ns
=========================================================  =========================

A value with units must have the right dimension (``--cutoff "5 ps"`` is
rejected). Values containing a space need quotes.

**In the Python API**, lengths and times can be pint Quantities, created with
the package registry:

.. code-block:: python

   from martinisoup.units import Q_
   from martinisoup.contact import ContactCounts

   cc = ContactCounts(proteins, metabolites, protein_types, cutoff_radius=Q_(0.5, "nm"))

A plain number passed to an analysis is taken to be in the **Universe's own
units**: Å (or ps) for a normally read Universe. The ``*_parallel`` and
``track_*`` wrappers read the files with default settings, so for them a
plain number is in Å.

Units in saved data
-------------------

Every pickle MartiniSoup writes holds plain numpy arrays and numbers, plus a
``units`` dict that maps saved keys to unit strings that pint can parse:

.. code-block:: python

   {
       'r': array([...]),
       'gr': array([...]),
       'mean_vol': 1.2e6,
       'units': {'r': 'angstrom', 'gr': 'dimensionless', 'mean_vol': 'angstrom ** 3', ...},
   }

The rules are:

- A unit given for a key applies to **every number in that value**: each
  element of an array, or each leaf of a nested dict.
- For a list of per-frame records whose fields have different units, the key
  is a **dotted path**: ``'frames.time'`` is the unit of the ``'time'`` field
  of every entry in ``data['frames']``.
- For the ``(values, counts)`` pairs in the ``residence-times --summary``
  output, ``'residences.values'`` and ``'residences.counts'`` give the unit of
  each half of the pair.
- Dimensionless numbers are listed explicitly as ``'dimensionless'``.

Reading saved data
^^^^^^^^^^^^^^^^^^

Without pint, read the unit strings directly. They use full unit names
(``angstrom``, ``picosecond``, ``centimeter ** 2 / second``).

With pint, rebuild Quantities in one line:

.. code-block:: python

   import pickle
   from martinisoup.units import get_quantity

   with open("protein_rdf.pkl", "rb") as fh:
       data = pickle.load(fh)

   r = get_quantity(data, "r")           # pint Quantity, e.g. in angstrom
   r_nm = r.to("nm").magnitude

or, in another program that has its own pint registry:

.. code-block:: python

   import pint
   ureg = pint.UnitRegistry()
   r = ureg.Quantity(data["r"], data["units"]["r"])

MartiniSoup uses a single pint registry, ``martinisoup.units.ureg``, which is
also set as pint's application registry. Use ``martinisoup.units.Q_`` (or
``ureg``) when combining your own Quantities with MartiniSoup's, because pint
does not allow Quantities from different registries to be mixed.

Files from older versions
^^^^^^^^^^^^^^^^^^^^^^^^^

Files written before units were recorded have no ``units`` entry. The fitters
still read them, with a warning, assuming the units that MartiniSoup has
always written: Å and ps (for ``msd-fitter``, any ``trajectory_units`` entry
in these files is ignored, because it describes the trajectory file, not the
data).

Units of the fitted results
---------------------------

- ``msd-fitter`` converts each input file from its recorded units to Å and ps
  before fitting, so replicas saved in different units are combined
  correctly. Diffusion coefficients are reported in cm²/s.
- ``residence-fitter`` converts lifetimes from their recorded units to ns.
  It refuses durations counted in frames (from ``use_time=False`` in the
  Python API), because these have no time unit to convert from.
