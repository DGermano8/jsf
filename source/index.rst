.. Jump-Switch-Flow documentation master file, created by
   sphinx-quickstart on Tue Dec  5 12:01:50 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the Jump-Switch-Flow Documentation!
==============================================

.. contents::
   :local:
   :depth: 2


Overview
--------

This package provides an algorithm for sampling from the
Jump-Switch-Flow (JSF) process. The JSF process is a continuous-time
process that can be used so represent compartmental models where
stochastic effects are important at low population sizes but can be
ignored at high population sizes.


Installation
------------

.. _installation:

This package is not yet available on PyPI. You can install it from a
local copy or from GitHub.

From Local Copy
^^^^^^^^^^^^^^^

If you have a local copy of the package, you can install it with pip.

.. code-block:: sh

   pip install /path/to/package

From GitHub
^^^^^^^^^^^

This won't work until the package has been made public. Once it has,
you can install it with pip.

.. code-block:: sh

   pip install git+https://github.com/DGermano8/jazz-shrill-fart.git

Testing
-------

There are some unit tests in the ``tests`` directory. You can run them
with the following command.

.. code-block:: sh

   python3 -m unittest discover -s tests

Housekeeping
------------

This package uses ``black`` and ``mypy`` for code formatting and type
checking, respectively. You can run them with the following commands.

.. code-block:: sh

   black jsf
   mypy jsf

..  LocalWords:  JSF
