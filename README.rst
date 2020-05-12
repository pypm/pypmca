====
pyPM.ca
====


.. image:: https://img.shields.io/pypi/v/pypmca.svg
        :target: https://pypi.python.org/pypi/pypmca

.. image:: https://img.shields.io/travis/deankarlen/pypmca.svg
        :target: https://travis-ci.com/pypmca/pypmca

* Free software: GNU General Public License v3


The ``pyPM.ca`` population modeller (www.pyPM.ca) describes connected systems with
discrete-time difference equations. It was developed specifically
to understand and characterize the CoViD-19 epidemic.

A ``pyPM.ca`` model is built by connecting a set of population objects with
connector objects. The connectors represent either a transfer that occurs
immediately at the next time step
or one which is delayed and distributed in time. Each
population object retains a record of its size at each time step, and also
maintains a list of future contributions, arising from delayed transfers
from other populations.
Calculations of population size are done either in terms of expectation values
or simulated data,
allowing the model to be used for both analysis and simulation.

Steady state solutions for these systems may develop. For viral epidemics these
consist of either exponential growth or decline of population size.
Systems perturbed by external influences, such as an instantaneous change
to a growth parameter or a sudden change in the size of a group,
may take time to relax to a new steady state solution.
The relaxation time depends on the various time delay distributions in the system.
To account for such perturbations, ``pyPM.ca`` models can include transition objects.

In order to model both long term and short term behaviour correctly,
realistic time delay distributions must be included.
The is achieved in ``pyPM.ca`` by using discrete-time difference equations,
in contrast to the ordinary differential equations (ODEs) approach
which form the foundation of the vast majority of epidemiological models.
The ODE approach is limited by its capability to introduce realistic time delays
and it offers no real benefit in modeling slowly evolving systems.

The ``pyPM.ca`` object oriented design makes it possible to create or modify a
connection network with little or no programming required, by using a suitable GUI.

A separate package, ``ipyPM``, provides a graphical user interface to ``pyPM.ca``
that runs inside a jupyter notebook. Data-model comparison, parameter adjustment,
and parameter estimation can be performed.
With these tools, the past behaviour of the system, and possible future
behaviours can be explored.
All aspects of the model can be adjusted, including connection network that
is at the core of the model.
This basic GUI can be used as an example for the design of special purpose GUIs
for other applications.

A single ``pyPM.ca`` model describes a homogenous system.
Heterogeneous systems can be modelled with an ensemble object,
by combining several ``pyPM.ca`` models, each
representing a distinct group or category, along with a contact matrix the represents
the engagement between the groups. Since an ensemble is a model object, tools that
interact with a model can also be used to interact with an ensemble.

Author
--------

 - `Dean Karlen <https://www.uvic.ca/science/physics/vispa/people/faculty/karlen.php>`_

Contributor
-------------

 - `Pradeep Reddy Raamana <https://crossinvalidation.com>`_



