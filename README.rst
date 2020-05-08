====
pypm
====


.. image:: https://img.shields.io/pypi/v/pypm.svg
        :target: https://pypi.python.org/pypi/pypm

.. image:: https://img.shields.io/travis/deankarlen/pypm.svg
        :target: https://travis-ci.com/pypm/pypm

* Free software: GNU General Public License v3


The ``pyPM`` population modeller (pyPM.ca) that describes connected systems with
discrete-time difference equations. It was developed specifically
to understand and characterize the CoViD-19 epidemic.

A ``pyPM`` model is built by connecting a set of population objects with
connector objects. The connectors represent either a transfer that is immediate
or one which is delayed and distributed in time. Each
population object retains a record of its number at each time step, and also
maintains a list of future contributions, arising from delayed transfers
from other populations.
Calculations are done either in terms of expectation values or simulated data,
allowing the model to be used as both an analysis tool and a simulator.

Steady state solutions for these systems may develop. For viral epidemics these
consist of either exponential growth or decline.
Systems perturbed by external influences, such as an instantaneous change
to a growth parameter or a sudden change in the size of a group,
may relax to a new steady state solution.
The relaxation time depends on the various time delay distributions in the system.
To account for such perturbations, ``pyPM`` models can include transition objects.

In order to correctly model both long term and short term behaviour, realistic time
delay distributions must be included.
The is achieved in ``pyPM`` by using discrete-time difference equations,
instead of ordinary differential equations (ODEs) which are the foundation of
the vast majority of epidemiological models.
The ODE approach is limited in its approach to introduce delays
and offer no real benefit in modeling slowing evolving systems.

The object oriented design makes it possible to create or modify a
connection network with little or no programming required.

A separate package, ``ipyPM``, provides a graphical user interface to pypm that runs
inside a jupyter notebook. Data-model comparison, parameter adjustment,
and parameter estimation can be performed.
With these tools, the past behaviour of the system, and possible future
behaviours can be explored.

A single ``pyPM`` model describes a homogenous system.
Heterogeneous systems can be modelled with an ensemble object,
by combining several ``pyPM`` models, each
representing a distinct group or category, along with a contact matrix the represents
the engagement between the groups.

Author
--------

 - `Dean Karlen <https://www.uvic.ca/science/physics/vispa/people/faculty/karlen.php>`_

Contributor
-------------

 - `Pradeep Reddy Raamana <https://crossinvalidation.com>`_



