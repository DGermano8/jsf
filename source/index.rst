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

JSF Mathematical Framework
--------------------------

To couple both the stochastic (Jumping) and deterministic (Flowing)
compartments, we model each compartment as to where they are in state space.

Consider a compartmental model with :math:`n` compartments :math:`\vec{V} = \left\{ V_i\right\}_{i=1}^n`, 
where :math:`V_{i}(t)` represents the value of the :math:`i` th compartment at time :math:`t`. 
For example, :math:`V_i` could be the number of people infected with a pathogen, or the copy 
number of a molecule in a cell. The state variables :math:`V_i` may take values from different 
domains depending upon the resolution needed for the model. For example, in an ODE, :math:`\vec{V}` 
will have real values and in a CTMC :math:`\vec{V}`  might have integer values.

Typically, discrete values are used to represent small populations, while larger populations 
will be represented with a continuum. To accommodate both scales, we model the domain of :math:`V_{i}`
as :math:`\mathcal{V}_{\Omega_i}=\{0,1,\ldots,\Omega_{i}\}\cup(\Omega_{i},\infty)`. The switching 
threshold parameter, :math:`\Omega_i\in \mathbb{Z}_{\geq 0}`, is where the :math:`i` th compartment 
transitions from discrete to continuous dynamics. If a compartment :math:`V_{i}` has a value in 
:math:`\{0,1,\ldots,\Omega_{i}\}`, we call it "discrete" (or "jumping"), and if it has a value in 
:math:`(\Omega_{i},\infty)`, we call it "continuous" (or "flowing"). While the switching threshold 
can be compartment specific, for ease of exposition, we will only consider a single threshold shared 
between all compartments :math:`\Omega = \Omega_i`. At any moment in time let us assume :math:`q` of the 
:math:`n` compartments are flowing. We use the notation :math:`\vec{V}_F = \left\{ V_i: V_i>\Omega \right\} \in (\Omega,\infty)^q`
and :math:`\vec{V}_J = \left\{ V_i: V_i\leq \Omega \right\} \in \left\{0, 1, \ldots, \Omega \right\}^{(n-q)}` 
to represent the compartments in each of the flowing and jumping states, respectively. 

The dynamics of each compartment :math:`V_i` are described by a set of :math:`m` "reactions" 
:math:`\mathcal{R} = \left\{\mathcal{R}_k \right\}_{k=1}^m`. Each reaction :math:`\mathcal{R}_k`
is defined by two properties: the rate (per unit time) at which it occurs, :math:`\lambda_{k}`, which may be 
(and usually is) a function of the state :math:`\vec{V}`; and the effect on the state, i.e. the change :math:`\eta_{ik}` 
to the size of compartment :math:`V_i` when reaction :math:`\mathcal{R}_k` occurs. As a matrix, 
:math:`\eta\in \mathbb{Z}^{n,m}` is referred to as the "stoichiometric matrix". For ODE models, these 
reactions occur continuously and are written in the form

.. math::
  \frac{\mathrm{d}\vec{V}}{\mathrm{d}t} = \eta \vec{\lambda}(\vec{V}),

while for CTMC models, reactions in the system :math:`\mathcal{R}` occur as discrete events.
In the later case, each reaction :math:`\mathcal{R}_k` has a separate propensity described by 
:math:`\lambda_k(\vec{V})`, this propensity remains constant between events but when an event 
:math:`\mathcal{R}_k` occurs, there is a change in :math:`\vec{V}` (as specified by the elements of 
:math:`\eta_{\cdot k}`) and therefore in :math:`\vec{\lambda}(\vec{V})`.

We define the subset of reactions :math:`\mathcal{S}\subseteq \mathcal{R}` to contain those treated 
as stochastic events. We define :math:`\mathcal{S}`, which we use throughout this manuscript, captures 
a larger set of reactions; :math:`\mathcal{S} = \left\{\mathcal{R}_k:\exists i \text{ s.t. } V_i\in\vec{V}_J \text{ and }  \left(\eta_{ik}\neq 0 \text{ or } \partial_{V_i}\lambda_k \neq 0\right)   \right\}`.
In this definition, a reaction is included in :math:`\mathcal{S}` if either (1) it causes a change in 
jumping (discrete) populations \textit{or} (2) it is influenced by a discrete population (perhaps as reactants for example). 


Reactions in :math:`\mathcal{S}` are simulated using stochastically sampled times similar to CTMC models. 
It is important to note that, unlike time homogeneous CTMC models, the propensities are not constant 
because the state :math:`\vec{V}` (and therefore :math:`\vec{\lambda}`) are continuously varying.
When any reaction :math:`\mathcal{R}_k\in\mathcal{S}` occurs, we say the system has "jumped" and an 
instantaneous change of :math:`\eta_{ik}` for each compartment :math:`V_i` occurs (irrespective of whether
:math:`V_i\in\vec{V}_J` or :math:`V_i\in\vec{V}_F` to ensure mass conservation is observed). We will refer 
therefore to reactions in :math:`\mathcal{S}` as "jumps". The reactions in :math:`\mathcal{S}'=\mathcal{R}\setminus \mathcal{S}`
are not stochastic, we call these "flows" because they represent the continual change of value of the 
relevant compartments, all of which are continuous by definition of :math:`\mathcal{S}'`.
At any moment in time, we denote :math:`|\mathcal{S}'| = p = m - |\mathcal{S}|` to be the number of reactions
which are flowing.

Finally, the hybrid model that we propose is capable of "switching". Switch events are defined as a 
compartment between :math:`\vec{V}_F` and :math:`\vec{V}_J`. These events occur when a compartment's 
value crosses the switching threshold :math:`\Omega`. Importantly, switch events can change 
:math:`\mathcal{S}` and are paradigm defining events which should occur infrequently compared to
jumps (frequent) and flows (continuous).

Due to the way that :math:`\mathcal{R}` is partitioned, it is possible to order the rows and columns of :math:`\eta` at 
any moment into the upper-triangular block form

.. math::
   \eta = \left(\begin{array}{c|c}
   \eta_{\mathcal{S}'} & \bar{\eta}_{\mathcal{S}} \\ \hline
   0 & \eta_{\mathcal{S}} \end{array}\right),

where :math:`\eta_{\mathcal{S}'} \in \mathbb{Z}^{q\times p}`, :math:`\eta_{\mathcal{S}} \in \mathbb{Z}^{(n-q)\times (m-p)}`
and :math:`\bar{\eta}_{\mathcal{S}} \in \mathbb{Z}^{q\times (m-p)}` refer to stoichiometric coefficients for changes in flowing 
compartments under flows, jumping compartments under jumps, and flowing compartments under jumps, respectively. Written 
as a system of equations analogous to (\ref{ODEmodels}), the hybrid JSF model we propose formally takes the following form. 
For any time interval :math:`t_0<t<t_1` between switching events, 

.. math::
    \frac{\mathrm{d} \vec{V}_F}{\mathrm{d} t} &= \eta_{\mathcal{S}'} \vec{\lambda}_{\mathcal{S}'}(\vec{V}) + \bar{\eta}_{\mathcal{S}} \vec{\Lambda}_{\mathcal{S}}(\vec{V}),\\
    \vec{V}_J(t) &= \vec{V}_J(t_0) + \eta_{\mathcal{S}} \int_{t_0}^t  \vec{\Lambda}_{\mathcal{S}}(\vec{V}) \ \mathrm{d} s,

where :math:`\vec{\lambda}_{\mathcal{S}'}\in\mathbb{R}^p` are the reaction rates of flows and :math:`\vec{\Lambda}_{\mathcal{S}}`
is a stochastic vector of :math:`m-p` delta-function spike trains that are derived from the realisations of :math:`m-p` 
different jumps sampled at rates which are dependent on the dynamic changes in the propensities :math:`\vec{\lambda}_{\mathcal{S}}\in\mathbb{R}^{m-p}` for these jumps 

The below figure shows how it is
possible for a variable to `switch` between flowing and jumping
regimes. When a flowing variable decreases to :math:`\Omega`, it switches
to jumping and we consider it a discrete variable. When a jumping
variable jumps from :math:`\Omega` to :math:`\Omega+1` it switches to flowing
and we consider it to be a continuous variable.

.. image:: _static/compartment_switching_depiction.png
   :width: 500
   :align: center
   :alt: compartment_switching_depiction

Lotka-Volterra preditor-prey model example
------------------------------------------

As a simple example for how JSF can be used to capture both stochastic and deterministic
dynamics, we consider the classic Lotka-Volterra preditor-prey model. 
In this model, individuals are either prey (:math:`V_1`) or preditors (:math:`V_2`).
This model is famous for being susceptible to Atto-fox problem, where the deterministic description
of the model allows for states to become infeasibly small, where the compartment would have otherwise had gone 
extinct. However, we will see that the JSF process can capture the stochastic effects of the model and permit 
the system to exhibit both the typical coexistence of the two species and the extinction of one of the species.

The model is described as the following way: the prey reproduce at a constant rate and are eaten by preditors at a rate proportional
to the number of preditors. The preditors die at a constant rate and are reproduced at 
a rate proportional to the number of prey they eat. Mathematically, the model is given by:

.. math::
   \frac{\mathrm{d} V_1}{\mathrm{d} t} &= \alpha V_1 - \beta V_1 V_2,\\
   \frac{\mathrm{d} V_2}{\mathrm{d} t} &=  \beta V_1 V_2 - \gamma V_2.


The compartment model diagram for this preditor-prey model is shown below.

.. image:: _static/preditorPreyModel.png
   :width: 400
   :align: center
   :alt: preditorPreyModel


Written as a JSF model, the stoichiometric matrices for the preditor-prey model is given by:

.. image:: _static/JSF_PredPrey.png
   :width: 500
   :align: center
   :alt: JSF_PredPrey


We can simulate this model using the JSF process. In doing so, we can see how the model behaves
as both a stochastic and deterministic process to obtain both the limit cycle and absorbing state behaviour, 
as shown in the figure below.

.. image:: _static/PP_Behaviour.png
   :width: 600
   :align: center
   :alt: PP_Behaviour


Epidemic simulation example
---------------------------

As a simple example, consider the SIS epidemic model. In this model, individuals are either susceptible 
(S) or infected (I). Susceptible individuals become infected at a rate proportional to the number of 
infected individuals, and infected individuals recover at a constant rate.

We can simulate this model using the JSF process. The only package that is needed is jsf, but we will 
import a few others to help us visualise the results.

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

.. Building the documentation
.. ^^^^^^^^^^^^^^^^^^^^^^^^^^

.. .. code-block:: sh

..    make html
..    cp build/html <my/website>

..  LocalWords:  JSF
