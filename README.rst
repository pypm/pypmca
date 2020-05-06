====
pypm
====


.. image:: https://img.shields.io/pypi/v/pypm.svg
        :target: https://pypi.python.org/pypi/pypm

.. image:: https://img.shields.io/travis/deankarlen/pypm.svg
        :target: https://travis-ci.com/pypm/pypm

* Free software: GNU General Public License v3


The ``pyPM`` population modeller (pyPM.ca) is a general framework for building models
of virus spread using discrete-time difference equations.


A ``pyPM`` model consists of a set of population objects connected by connector
objects, allowing expectations and simulated data for time histories (such as
positive tests and ICU admissions) to be produced. Propagation from one
population to another is done with realistic time lags, defined by normal delay
distributions, with the mean and standard deviation specified for each connector.
It is important to note that models based on ODEs (ordinary differential
equations), popular in the virus modelling community perhaps for historical
reasons, are not able to correctly model the propagation of infectious bursts to
the reporting stage.


Authors
--------

 - `Dean Karlen <https://www.uvic.ca/science/physics/vispa/people/faculty/karlen.php>`_

Contributors
-------------

 - `Pradeep Reddy Raamana <https://crossinvalidation.com>`_



