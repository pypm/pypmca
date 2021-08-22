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

* 0.2.7 (2020-11-21)

 * condition simulation ensembles

* 0.2.8 (2020-11-28)

 * add tool to make intervals

* 0.2.9 (2020-11-29)

 * minor corrections and improvements

* 0.2.10 (2020-12-04)

 * allow gamma delay distributions
 * new reference model
 * smearing parameters for intervals
 * weekly noise reporting

* 0.2.11 (2020-12-12)

 * produce intervals for multiple categories at once

* 0.2.12 (2020-12-27)

 * add AgeModeller tool to automatically fit age groups
 * improve functionality of IntervalMaker
 * fix chi^2 for cumul_reset = True

* 0.2.13 (2020-12-29)

 * modify linear modifier
 * fix chi^2 for daily data
 * new reference model

* 0.2.15 (2020-12-30)

 * minor corrections and additions

* 0.2.16 (2021-01-02)

 * fix for interval maker for scenarios with future transitions
 * modified file format for scenario hub

* 0.2.17 (2021-01-27)

 * fixes for backup parameter for interval maker
 * set protocol level for pickling models
 * add new model (2_8) with second variant infection cycle

* 0.2.18 (2021-01-31)

 * let fit skip specified dates (xmas period for example)

* 0.2.19 (2021-02-14)

 * add option to produce simulated weekly data
 * minor fix for copy

* 0.2.20 (2021-03-10)

 * fix issues with Ensemble

* 0.2.21 (2021-03-28)

 * fix mistake in Ensemble transitions

* 0.2.22 (2021-08-22)

 * add feature to combine evolve_expectations and generate_data
 * improve optimization for local fits (cumul_reset)