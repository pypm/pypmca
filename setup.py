#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

descr = """
The pyPM population modeller (pyPM.ca) is a general framework for building models 
of virus spread using discrete-time difference equations.  

A pyPM model consists of a set ofpopulation objects connected by connector 
objects, allowing expectations and simulateddata for time histories (such as 
positive tests and ICU admissions) to be produced. Propa-gation from one 
population to another is done with realistic time lags, defined by normaldelay 
distributions, with the mean and standard deviation specified for each connector. 
Itis important to note that models based on ODEs (ordinary differential 
equations), popularin the virus modelling community perhaps for historical 
reasons, are not able to correctlymodel the propagation of infectious bursts to 
the reporting stage.
"""

requirements = [ 'scipy', 'numpy' ]

setup_requirements = ['pytest-runner', ] + requirements

test_requirements = ['pytest>=3', ] + requirements

setup(
    author="Dean Karlen",
    author_email='karlen@uvic.ca',
    python_requires='>=3.5',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="population modeller",
    install_requires=requirements,
    license="GNU General Public License v3",
    long_description=descr,
    include_package_data=True,
    keywords='pypm',
    name='pypm',
    packages=find_packages(include=['pypm', 'pypm.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/pypm/pypm',
    version='0.0.5',
    zip_safe=False,
)
