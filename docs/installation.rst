Installation
============

Requirements
------------

- Python >= 3.11
- `uv <https://docs.astral.sh/uv/>`_ (see the `uv installation guide
  <https://docs.astral.sh/uv/getting-started/installation/>`_)
- Dependencies (installed automatically): numpy, MDAnalysis >= 2.8, lmfit,
  scipy, freud-analysis, pandas, requests, tqdm, tidynamics

Install
-------

From the repository root:

.. code-block:: bash

   uv sync --no-dev

This creates a virtual environment in ``.venv`` and installs MartiniSoup with
the exact dependency versions recorded in ``uv.lock``, so every install of a
given commit gets the same environment. Either activate it
(``source .venv/bin/activate``) or prefix commands with ``uv run``, e.g.
``uv run martinisoup --help``.

Development
-----------

.. code-block:: bash

   uv sync

This installs MartiniSoup in editable mode (changes to ``src/`` take effect
without reinstalling) together with the development dependencies, such as
pytest. Run the tests with:

.. code-block:: bash

   uv run pytest

To add or upgrade a dependency, use ``uv add <package>`` (or
``uv add --dev <package>`` for development-only tools) and commit the
updated ``pyproject.toml`` and ``uv.lock`` together.

Command-line only
-----------------

If you only need the ``martinisoup`` command and not the Python API, install
it as a standalone tool in its own isolated environment:

.. code-block:: bash

   uv tool install .

or directly from GitHub:

.. code-block:: bash

   uv tool install git+https://github.com/Martini-Force-Field-Initiative/MartiniSoup.git

``uv tool install`` resolves dependencies from ``pyproject.toml`` and does
not use ``uv.lock``.
