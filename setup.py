# -*- coding: utf-8 -*-
# Copyright (C) 2019  Thibault Hallouin
from setuptools import setup


with open("README.md", "r") as fh:
    long_desc = fh.read()

with open('eflowcalc/version.py') as fv:
    exec(fv.read())

setup(
    name='eflowcalc',

    version=__version__,

    description='EFlowCalc: An Efficient Calculator of Ecological Streamflow Characteristics in Python',
    long_description=long_desc,
    long_description_content_type="text/markdown",

    url='https://github.com/ThibHlln/eflowcalc',

    author='Thibault Hallouin',
    author_email='thibault.hallouin@ucdconnect.ie',

    license='GPLv3',

    classifiers=[
        'Development Status :: 4 - Beta',

        'Natural Language :: English',

        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Hydrology',

        'Operating System :: MacOS',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',


        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: Implementation :: CPython'
    ],

    packages=['eflowcalc'],

    install_requires=[
        'numpy',
        'pandas'
    ],
)
