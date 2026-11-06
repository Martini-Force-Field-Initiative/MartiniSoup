"""
Unit handling for martinisoup, built on pint.

The package uses one shared pint registry, ``ureg``, which is also set as
pint's application registry. A single registry matters: pint cannot add,
compare or unpickle together Quantities that come from different registries.

Units of a Universe
-------------------
MDAnalysis converts coordinates, box dimensions and times to its own base
units (Å, ps) as it reads a file, unless the reader was created with
``convert_units=False``. ``u.trajectory.units`` describes the *file* (e.g.
``nm`` for .xtc), not the data in memory, so it must not be used on its own.
:meth:`UnitSystem.from_universe` applies the correct rule.

Saved data convention
---------------------
Programs save plain numpy arrays and numbers, plus a ``'units'`` dict mapping
saved keys to pint-parseable unit strings (e.g. ``'angstrom ** 2'``):

- A unit given for a key applies to every number in that value: each
  element of an array, or each leaf of a nested dict (e.g. all counts in
  ``{protein: {metabolite: {'count': ..., 'normalised_count': ...}}}``).
- For a list of per-frame records (dicts) whose fields have different
  units, the key is a dotted path, e.g. ``'frames.time'`` for the ``'time'``
  field of every entry in ``data['frames']``.
- Dimensionless numbers are listed explicitly as ``'dimensionless'``.

Consumers without pint can read the strings; pint users can rebuild
Quantities with :func:`get_quantity`, or directly with
``ureg.Quantity(data[key], data['units'][key])``.

Numbers without units
---------------------
- Python API: a plain number passed where a length or time is expected is
  taken to be in the *Universe's* units (Å / ps for a Universe read with
  default settings). Pass a pint Quantity to use any other unit.
- Command line: a bare number is in the unit stated in each option's help
  (Å for distances); a value such as ``"0.5 nm"`` is converted.
"""

from __future__ import annotations

import argparse
import warnings
from dataclasses import dataclass

import numpy as np
import pint
from MDAnalysis import units as mda_units

ureg = pint.UnitRegistry()
# CHARMM's AKMA time unit, which MDAnalysis can report for some trajectories.
# Written out (rather than read from MDAnalysis.units at import) so this module
# imports even where MDAnalysis is mocked, e.g. the docs build; the value is
# checked against MDAnalysis's factor in tests/test_units.py.
ureg.define("AKMA_time = 4.888821e-2 * picosecond")
pint.set_application_registry(ureg)
Q_ = ureg.Quantity

UNITS_KEY = 'units'

# MDAnalysis unit names -> pint unit names. Explicit because some MDAnalysis
# names mean something else to pint (e.g. 'A' is ampere in pint).
_MDA_TO_PINT = {
    # length (MDAnalysis.units.lengthUnit_factor)
    'A': 'angstrom', 'Angstrom': 'angstrom', 'angstrom': 'angstrom',
    'Å': 'angstrom',  # ANGSTROM SIGN, the character MDAnalysis uses
    'Å': 'angstrom',  # LATIN CAPITAL LETTER A WITH RING ABOVE
    'nm': 'nanometer', 'nanometer': 'nanometer',
    'pm': 'picometer', 'picometer': 'picometer',
    'fm': 'femtometer', 'femtometer': 'femtometer',
    # time (MDAnalysis.units.timeUnit_factor)
    'ps': 'picosecond', 'picosecond': 'picosecond',
    'fs': 'femtosecond', 'femtosecond': 'femtosecond',
    'ns': 'nanosecond', 'nanosecond': 'nanosecond',
    'us': 'microsecond', 'μs': 'microsecond', 'microsecond': 'microsecond',
    'ms': 'millisecond', 'millisecond': 'millisecond',
    's': 'second', 'sec': 'second', 'second': 'second',
    'AKMA': 'AKMA_time',
}


def pint_unit_name(mda_name: str) -> str:
    """Translate an MDAnalysis unit name (e.g. ``'A'``, ``'nm'``) to pint's."""
    try:
        return _MDA_TO_PINT[mda_name]
    except KeyError:
        raise ValueError(f"Unknown MDAnalysis unit {mda_name!r}") from None


def unit_string(unit) -> str:
    """Canonical, pint-parseable string for a unit, Quantity or unit string."""
    if isinstance(unit, pint.Quantity):
        unit = unit.units
    return str(ureg.Unit(unit) if isinstance(unit, str) else unit)


def pretty(unit) -> str:
    """Short display form of a unit for plot labels and CSV headers, e.g. 'Å²'."""
    u = ureg.Unit(unit) if isinstance(unit, str) else unit
    return f"{u:~P}" if u != ureg.dimensionless else ""


@dataclass(frozen=True)
class UnitSystem:
    """
    Length and time units of the data in a Universe.

    Stored as unit strings so the object pickles cleanly to worker processes.
    """

    length: str = 'angstrom'
    time: str = 'picosecond'

    @classmethod
    def from_universe(cls, universe) -> "UnitSystem":
        """
        Units of the coordinates, box and times held in memory by ``universe``.

        If the reader converted on load (the MDAnalysis default) the data is
        in MDAnalysis base units (Å, ps). Only if it was created with
        ``convert_units=False`` are ``trajectory.units`` (the file's native
        units) the units of the data. Readers with no conversion flag or no
        recorded units (e.g. in-memory trajectories) hold base units.
        """
        return cls.from_trajectory(universe.trajectory)

    @classmethod
    def from_trajectory(cls, reader) -> "UnitSystem":
        """As :meth:`from_universe`, given the trajectory reader directly."""
        base = mda_units.MDANALYSIS_BASE_UNITS
        length, time = base['length'], base['time']
        if not getattr(reader, 'convert_units', True):
            native = getattr(reader, 'units', {}) or {}
            length = native.get('length') or length
            time = native.get('time') or time
        return cls(length=pint_unit_name(length), time=pint_unit_name(time))

    @property
    def length_unit(self) -> pint.Unit:
        return ureg.Unit(self.length)

    @property
    def time_unit(self) -> pint.Unit:
        return ureg.Unit(self.time)

    def length_magnitude(self, value):
        """
        Magnitude of a length in this system's length unit.

        A pint Quantity is converted; a plain number is taken to already be
        in this system's length unit.
        """
        return _magnitude(value, self.length_unit, 'length')

    def time_magnitude(self, value):
        """Magnitude of a time in this system's time unit (see length_magnitude)."""
        return _magnitude(value, self.time_unit, 'time')


def _magnitude(value, unit, kind):
    if isinstance(value, pint.Quantity):
        if not value.check(unit.dimensionality):
            raise pint.DimensionalityError(value.units, unit, extra_msg=f" (expected a {kind})")
        return value.to(unit).magnitude
    return value


# --------------------------------------------------------------------------- #
# Command-line parsing
# --------------------------------------------------------------------------- #

def quantity_arg(default_unit: str):
    """
    argparse ``type=`` for an option that takes a value with optional units.

    A bare number (``"5"``) is taken to be in ``default_unit``; a value with
    units (``"0.5 nm"``) must have the same dimensionality and is kept as
    given. Returns a pint Quantity.
    """
    target = ureg.Unit(default_unit)

    def parse(text: str) -> pint.Quantity:
        try:
            q = Q_(text)
        except Exception as exc:
            raise argparse.ArgumentTypeError(f"cannot parse {text!r} as a quantity: {exc}")
        if q.unitless:
            return Q_(q.magnitude, target)
        if not q.check(target.dimensionality):
            raise argparse.ArgumentTypeError(
                f"{text!r} has units of {q.units:~P}, expected a quantity convertible to {target:~P}"
            )
        return q

    parse.__name__ = f"quantity in {target}"
    return parse


# --------------------------------------------------------------------------- #
# Reading saved data
# --------------------------------------------------------------------------- #

def get_units(data: dict, legacy: dict | None = None, source: str = 'data') -> dict:
    """
    The ``'units'`` dict of saved data.

    Files written before units were recorded have no ``'units'`` key; for
    those, ``legacy`` gives the units the data is known to be in, and a
    warning is issued.
    """
    if UNITS_KEY in data:
        return data[UNITS_KEY]
    if legacy is None:
        raise KeyError(f"{source} has no {UNITS_KEY!r} entry")
    warnings.warn(
        f"{source} has no recorded units (written by an older martinisoup); "
        f"assuming {legacy}",
        stacklevel=2,
    )
    return legacy


def get_quantity(data: dict, key: str, units: dict | None = None):
    """Rebuild ``data[key]`` as a pint Quantity using the saved units."""
    units = units if units is not None else data[UNITS_KEY]
    return Q_(np.asarray(data[key]), units[key])
