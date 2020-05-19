#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages
import versioneer

descr = """
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
"""

requirements = ['scipy', 'numpy']

setup_requirements = ['pytest-runner', ] + requirements

test_requirements = ['pytest>=3', ] + requirements

setup(
        author="Dean Karlen",
        author_email='karlen@uvic.ca',
        python_requires='>=3.5',
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: End Users/Desktop',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Natural Language :: English',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            ],
        description="pyPM.ca population modeller",
        install_requires=requirements,
        license="GNU General Public License v3",
        long_description=descr,
        include_package_data=True,
        keywords='pypmca',
        name='pypmca',
        version=versioneer.get_version(),
        packages=find_packages(include=['pypmca', 'pypmca.*']),
        setup_requires=setup_requirements,
        test_suite='tests',
        tests_require=test_requirements,
        url='http://www.pypm.ca',
        zip_safe=False,
        cmdclass=versioneer.get_cmdclass()
        )
