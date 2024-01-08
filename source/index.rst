.. Jump-Switch-Flow documentation master file, created by
   sphinx-quickstart on Tue Dec  5 12:01:50 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the Jump-Switch-Flow Documentation!
==============================================

.. image:: _static/jsf-logo.png
   :width: 400
   :align: center

.. contents::
   :local:
   :depth: 1

Overview
--------

This package provides an algorithm for sampling from the
Jump-Switch-Flow (JSF) process. The JSF process is a continuous-time
process that can be used to represent compartmental models where
stochastic effects are important at low population sizes but can be
ignored at high population sizes.

Method
------

To couple both the stochastic (Jumping) and deterministic (Flowing)
compartments, we model each compartment as to where they are in state space.

Consider a compartmental model where :math:`V_{i}(t)` represents the value
of the :math:`i` th compartment at time :math:`t`. In a differential equation based
model, the :math:`V_{i}` could take values in :math:`\mathbb{R}_{\geq 0}`, and
in a CTMC model they might take values in :math:`\mathbb{Z}_{\geq 0}`.

To couple both the stochastic (Jumping) and deterministic (Flowing)
compartments in Jump-Switch-Flow, we model each compartment as taking
values in :math:`\mathcal{V}_{\Omega} = \{0,1,\ldots,\Omega\} \cup (\Omega,\infty)`.
The value :math:`\Omega` is the value at which the system transitions from discrete
to continuous dynamics. If a compartment :math:`V_{i}` has a value in
:math:`\{0,1,\ldots,\Omega\}`, we describe that compartment as
`discrete` (or `Jumping`), and if it has a value in
:math:`(\Omega,\infty)`, we describe it as `continuous` (or
`Flowing`).

Changes to the :math:`V_i` in these models are defined by a set of
`reactions`. Each reaction consists of three things: the rate it
occurs, the reactants consumed it, and the products generated. For a
given reaction, we refer to the difference in the amount produced, and
the mount consumed as the `flow`.

The set :math:`\mathcal{S}` consists of the reactions which have at least
one discrete reactant or product. We refer to the occurrence of these
reactions as a `jumps` because they involve one of the discrete
variables values changing discretely. The reactions not in
:math:`\mathcal{S}` are called `flows`` because they represent the
continual change of value.

Jump events occur following an inhomogeneous arrival process. Reaction
:math:`j` occurs at a rate :math:`\lambda_{j}`, which may depend upon the
values of all the reactants of that reaction. The net rate of
reactions is :math:`\sum_{j \in \mathcal{S}} \lambda_{j}`.
For a detailed description of how to sample these reaction times, see the description
of time-varying Poisson arrival processes by `Klein and Roberts <https://journals.sagepub.com/doi/10.1177/003754978404300406>`_.

For each possible reaction in the system, :math:`k`, we include
a new variable, :math:`J_{k}(t)`, which counts down the time until that
reaction will occur. At time :math:`t_0` the value of :math:`J_{k}(t_0)` is
initialised with a uniform random variable,
:math:`u_{k}\sim\text{Uniform(0,1)}`. The remaining time until the
reaction fires then decreases as described by the following equation:

.. math::
  J_{k}(t)=u_{k} - \left( 1 - \exp\left\{- \int_{t_{0}}^{t} \lambda_{k} ds \right\} \right),

where :math:`\lambda_{k}` is the rate of the :math:`k` th reaction type. Since
the reaction rates are continuous, this is a strictly decreasing
function of time. Once one of the :math:`J_k` reaches zero, that reaction
occurs and all of the :math:`J_i` are rest. Note that the rate of
reactions can change continuously through time. This is because the
rate may depend on other variables that are governed by differential
equations.

Between jump events, the discrete variables remain constant and the
change of continuous variables follows a system of differential
equations: :math:`dV_i/dt = \sum_{j\in \mathcal{S}^{c}} \lambda_j \Delta_{i,j}`,
where :math:`\Delta_{i,j}` is the change in the amount of
variable :math:`i` during a reaction :math:`j`.

The below figure shows how it is
possible for a variable to `switch` between flowing and jumping
regimes. When a flowing variable decreases to :math:`\Omega`, it switches
to jumping and we consider it a discrete variable. When a jumping
variable jumps from :math:`\Omega` to :math:`\Omega+1` it switches to flowing
and we consider it to be a continuous variable.

.. image:: _static/compartment_switching_depiction.png
   :width: 600
   :align: center
   :alt: SIS epidemic example

Epidemic simulation example
---------------------------

As a simple example, consider the SIS epidemic model. In this model,
individuals are either susceptible (S) or infected (I). Susceptible
individuals become infected at a rate proportional to the number of
infected individuals, and infected individuals recover at a constant
rate.

We can simulate this model using the JSF process. The only package
that is needed is `jsf`, but we will import a few others to help us
visualise the results.

.. code-block:: python

   import pandas as pd
   from plotnine import *
   import random
   import jsf
   random.seed(7)

Defining the SIS model
^^^^^^^^^^^^^^^^^^^^^^

Next, we define the initial condition of the SIS model and the
infection and recovery process. The SIS model has two compartments,
so we define the initial condition as a list of length two. The
infection and recovery process is defined by a function that takes
the current state of the system and the current time and returns a
list of length two containing the rates of infection and recovery.

The reactant and product matrices are used to define the stoichiometry
of the process. The reactant matrix defines the change in the number
of individuals in each compartment when a reaction occurs, and the
product matrix defines the change in the number of individuals in each
compartment after a reaction occurs. For the SIS model, the reactant
matrix is ``[[1, 1], [0, 1]]`` and the product matrix is
``[[0, 2], [1, 0]]``. This means that when an infection occurs, the
number of susceptible individuals decreases by one and the number of
infected individuals increases by one. When a recovery occurs, the
number of infected individuals decreases by one and the number of
susceptible individuals increases by one.

.. code-block:: python

   x0 = [1000 - 2, 2]
   rates = lambda x, _: [2e-3 * x[0] * x[1], 1.0 * x[1]]
   reactant_matrix = [[1, 1], [0, 1]]
   product_matrix = [[0, 2], [1, 0]]

Finally, we define the maximum time of the simulation.

.. code-block:: python

   t_max = 10.0

There is a little bit of configuration needed to tell JSF how to
actually run the simulation.

.. code-block:: python

   stoich = {
       "nu": [
           [a - b for a, b in zip(r1, r2)]
           for r1, r2 in zip(product_matrix, reactant_matrix)
       ],
       "DoDisc": [1, 1],
       "nuReactant": reactant_matrix,
       "nuProduct": product_matrix,
   }
   my_opts = {"EnforceDo": [0, 0], "dt": 0.1, "SwitchingThreshold": [50, 50]}

Then we can call `jsf.jsf` to simulate the process using the operator
splitting method.

.. code-block:: python

   sim = jsf.jsf(x0, rates, stoich, t_max, config=my_opts, method="operator-splitting")

Visualising the simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, we can plot the results of the simulation. We'll use a
combination of `pandas` and `plotnine` to do this, but the output of
`jsf` is a list of numbers so it should be easy to use whichever
plotting library you prefer.

.. code-block:: python

   sim_df = pd.DataFrame(
       {"time": sim[1], "susceptible": sim[0][0], "infectious": sim[0][1]}
   ).melt(id_vars=["time"], value_vars=["susceptible", "infectious"])

   sim_p9 = (
       ggplot()
       + geom_line(data=sim_df, mapping=aes(x="time", y="value", colour="variable"))
       + geom_hline(yintercept=my_opts["SwitchingThreshold"][1], linetype="dashed")
       + scale_y_sqrt(name="Population size")
       + labs(x="Time", colour="Status")
       + theme(legend_position="top")
       + theme_bw()
   )

   sim_p9.save("sis_example.png", height=4, width=6)

Which gives us the following plot. Note that initially the process is
stochastic as it jumps around before hitting the threshold at which
point it follows the differential equations.

.. image:: _static/sis_example.png
   :width: 700
   :align: center
   :alt: SIS epidemic example

Types
-----

The ``jsf.types`` module provides some key types for this package.
There is nothing fancy here; they are just used to make the type hints
more informative and help to leverage ``mypy``.

- ``CompartmentValue``: the value of a compartment, this is a ``float``.
- ``SystemState``: the state of the system, this is a list of ``CompartmentValue`` s.
- ``Time``: the time, this is a float.

Recall you can use the following to type check the code:

.. code-block:: sh

    mypy jsf tests


API
---

.. toctree::
   :maxdepth: 2
   :caption: Modules

   jsf
   exact

Installation
------------

.. _installation:

This package is not yet available on PyPI. You can install it from a :ref:`local copy <local_copy_installation>` or from :ref:`GitHub <github_installation>`.

.. _local_copy_installation:

From Local Copy
^^^^^^^^^^^^^^^

If you have a local copy of the package, you can install it with pip.

.. code-block:: sh

   pip install /path/to/package

.. _github_installation:

From GitHub
^^^^^^^^^^^

This won't work until the package has been made public. Once it has,
you can install it with pip.

.. code-block:: sh

   pip install git+https://github.com/DGermano8/jsf.git

FAQs
----

If you have a question that is not answered by this documentation,
please lodge an issue on the GitHub page for this package:
https://github.com/DGermano8/jsf

Housekeeping
------------

Testing
^^^^^^^

There are some unit tests in the ``tests`` directory. You can run them
with the following command.

.. code-block:: sh

   python3 -m unittest discover -s tests

Code formating and checking
^^^^^^^^^^^^^^^^^^^^^^^^^^^

This package uses ``black`` and ``mypy`` for code formatting and type
checking, respectively. You can run them with the following commands.

.. code-block:: sh

   black jsf
   mypy jsf

Building the documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: sh

   make html
   cp build/html <my/website>

..  LocalWords:  JSF
