====
pypm
====


.. image:: https://img.shields.io/pypi/v/pypm.svg
        :target: https://pypi.python.org/pypi/pypm

.. image:: https://img.shields.io/travis/deankarlen/pypm.svg
        :target: https://travis-ci.com/pypm/pypm

* Free software: GNU General Public License v3


The ``pyPM`` population modeller (pyPM.ca) is a general framework for population
modeling using discrete-time difference equations. It was developed specifically
to understand and characterize the CoViD-19 epidemic.

A ``pyPM`` model is built by connecting a set of population objects with
connector objects. The connectors represent either a transfer that is immediate
or one which is delayed and distributed in time. Each
population object retains a record of its number at each step, and also
a list of future contributions, arising from the delayed transfers from other
populations.
The calculations can be done in terms of expectation values or simulated data.

By using discrete-time difference equations, realistic delays are easily
incorporated.
This is in stark contrast with the majority of
epidemiological models, which are almost exclusively based on
ordinary differential equations (ODEs).
With a single stage between two populations, ODE delays follow an
exponential distribution.

The object oriented design makes it possible to model complex systems with little
or no programming required.

A separate package, ipypm, provides a graphical user interface to pypm that runs
inside jupyter notebook. Data-model comparison and parameter adjustment
can be performed.

Heterogeneous systems can be modelled by combining several models, each
representing a distinct group or category, into an ensemble.


Author
--------

 - `Dean Karlen <https://www.uvic.ca/science/physics/vispa/people/faculty/karlen.php>`_

Contributor
-------------

 - `Pradeep Reddy Raamana <https://crossinvalidation.com>`_



