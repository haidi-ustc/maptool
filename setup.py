#!/usr/bin/env python
# coding: utf-8

import sys
import platform

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext as _build_ext
default_prefix='maptool'

setup(
    name="maptool",
    packages=find_packages(),
    version="0.1.0",
    setup_requires=['pymatgen','dpdata'],
    install_requires=[],
    package_data={},
    author="haidi",
    author_email="haidi@mail.ustc.edu.cn",
    maintainer="haidi",
    maintainer_email="haidi@mail.ustc.edu.cn",
    url="https://github.com/haidi-ustc/maptools",
    license="MIT",
    description=default_prefix+" is a pre- and post-processing software for materials science"
                "This software is based on Pymatgen and some open-source Python library.",
    keywords=["VASP", "Lammps", "QuantumEspresso","ASE","Pymatgen",
              "electronic", "structure", "analysis", "phase", "diagrams"],
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules"
    ],
    entry_points={
          'console_scripts': [
              'mpt = maptool.mpt:main'],
 }
)
