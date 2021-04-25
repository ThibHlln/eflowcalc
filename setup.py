# Copyright (C) 2019-2021  Thibault Hallouin
from setuptools import setup


with open("README.rst", "r") as fh:
    long_desc = fh.read()

with open('eflowcalc/version.py') as fv:
    exec(fv.read())

setup(
    name='eflowcalc',

    version=__version__,

    description='EFlowCalc: A Calculator of Streamflow Characteristics in Python',
    long_description=long_desc,
    long_description_content_type="text/x-rst",

    download_url="https://pypi.python.org/pypi/eflowcalc",
    project_urls={
        "Bug Tracker": "https://github.com/thibhlln/eflowcalc/issues",
        "Source Code": "https://github.com/thibhlln/eflowcalc",
        "Documentation": "https://thibhlln.github.io/eflowcalc",
    },

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

        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: Implementation :: CPython'
    ],

    packages=['eflowcalc'],

    install_requires=[
        'numpy',
        'pandas'
    ],

    extras_require={
        'tests': ['netCDF4', 'cftime']
    }
)
