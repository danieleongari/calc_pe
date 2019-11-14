#!/usr/bin/env python
from __future__ import absolute_import
from setuptools import setup, find_packages

if __name__ == '__main__':
    setup(
        packages=find_packages(),
        long_description=open('README.md').read(),
        long_description_content_type='text/markdown',
        name='calc_pe',
        author='Daniele Ongari',
        author_email='daniele.ongari@epfl.ch',
        description=
        'Calculator for the CO2 parasitic energy from CO2 and N2 isotherms',
        url='https://github.com/danieleongari/calc_pe',
        license='MIT',
        classifiers=['Programming Language :: Python'],
        version='1.0.1',
        install_requires=['numpy', 'pandas>=0.24.0', 'pyiast'],
        scripts=['bin/calc_pe'],
    )
