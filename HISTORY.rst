=======
History
=======

0.1 (2020-05-12)
------------------

* First release on PyPI.
* 0.1.1 - minor fix to __init__.py
* 0.1.2 - bug fix in evolve_expectations
* 0.1.3 - remove min/max parameter check
* 0.1.4 - revised the likelihood calculation for MCMC (now shape and normalization treated separately)

0.2 (2020-05-19)
----------------

* 0.2.1 beta release on PyPI
  (includes new reference model 2.1, and copy methods)
* 0.2.2 (2020-06-16)

  * fix issue with in_icu for reference model #1
  * allow no reports on some days of week (eg. BC: Sundays)
  * add fit statistics to help tune reporting noise
  * fix bugs in Ensemble

* 0.2.3 (2020-07-30)

  * allow cumulative to start from 0 at first point in fit
  * allow modifier to change parameter linearly with time
  * add code for submitting .csv file for US Forecast Hub
  * do interpolation with log(alpha)
  * fix bugs

* 0.2.4 (2020-07-31)

  * minor bug fix

* 0.2.5 (2020-09-11)

  * bug fixes for ensemble
  * updates for forecast hub

* 0.2.6 (2020-10-04)

 * addition functionality to model vaccination