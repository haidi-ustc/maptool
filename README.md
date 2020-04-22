## Introduction
[![Build Status](https://travis-ci.org/haidi-ustc/maptool.svg?branch=master)](https://travis-ci.org/haidi-ustc/maptool)

Maptool (Material Project Tool) is an open source python package for pre- and post- process input and output of first principles simulation software, which is based on Pymatgen and ASE and other open source python packages. Currently, the VASP code are mainly supported, later it will support more DFT codes. The package allows for generating, building, modifying and analysing of both molecular and crystal structure. Besides, maptool can automatically generating the input files for different DFT code and analysing the output, all results are saved with pretty format and figure based on matplotlib. See the [documentation](https://github.com/haidi-ustc/maptool/tree/master/doc/) for information about installation and usage.

## Features
* Random generation of atomic (3D, 2D, and 1D) crystals and point group clusters for a given symmetry group and stoichiometry
* Building of cluster, shell strcuture and surface
* Add strain to crystal structure
* Analysing of symmetry of molecule and crystal structure
* Structure transformation to cif, xsf, xyz, mol or POSCAR files via pymatgen and ASE
* General operation for 2D materials (splitting, strain, bending and so on)
* VASP input generation and output analysing and visualization
* Online retrieval function based on pymatgen
* Local database for storing the previous result

## Dependencies
* [SciPy 1.3.1](https://www.scipy.org/install.html)
* [NumPy 1.17.3](https://www.scipy.org/scipylib/download.html)
* [Pandas 0.25.2](https://pandas.pydata.org/getpandas.html)
* [Pymatgen 2019.10.16](http://pymatgen.org/#getting-pymatgen)
* [ASE 3.19.0](https://pypi.org/project/ase)

### Optional
* [openbabel 2.4.1 (Python bindings)](http://openbabel.org/wiki/Main_Page) (allows for additional molecule file formats)
* [pymongo](https://api.mongodb.com/python) (allows for database operation)

## Installation

```
$ git clone git@github.com:haidi-ustc/maptool.git
$ cd maptool
$ python setup.py install
```
To use online retrieval function, the materials project `API_KEK` should be added into the ~/.bashrc. The `API_KEY` can be obtained from [materials project](https://www.materialsproject.org).  

```
export MAPI_KEY='fdafdwe203213faf22'
```
use your own key replace the above string and write it into ~/.bashrc file.

After installation, excute the following command.
```
mpt -i
```
it will shows similar information:
```
maptool
------------

Version: 0.1.dev47+g08014e6.d20200419
Path:    /home/dgx/software/maptools/maptool
Date:    Apr-19-2020

Dependency
------------
Python version=3.7.4 (default, Aug 13 2019, 20:35:49) 
[GCC 7.3.0]

   pymongo      3.9.0   /home/dgx/miniconda3/lib/python3.7/site-packages/pymongo
     numpy     1.17.3   /home/dgx/miniconda3/lib/python3.7/site-packages/numpy
     scipy      1.3.1   /home/dgx/miniconda3/lib/python3.7/site-packages/scipy
    mayavi      4.7.1   /home/dgx/miniconda3/lib/python3.7/site-packages/mayavi
matplotlib      3.1.1   /home/dgx/miniconda3/lib/python3.7/site-packages/matplotlib
      tqdm     4.36.1   /home/dgx/miniconda3/lib/python3.7/site-packages/tqdm
    dpdata     0.1.2    /home/dgx/miniconda3/lib/python3.7/site-packages/dpdata/dpdata
      nose      1.3.7   /home/dgx/miniconda3/lib/python3.7/site-packages/nose
  coverage        5.1   /home/dgx/miniconda3/lib/python3.7/site-packages/coverage
    spglib     1.14.1   /home/dgx/miniconda3/lib/python3.7/site-packages/spglib
    pyhull   2015.2.1   /home/dgx/miniconda3/lib/python3.7/site-packages/pyhull
  pymatgen 2019.10.16   /home/dgx/miniconda3/lib/python3.7/site-packages/pymatgen
      qmpy            Not Found
       ase     3.19.0   /home/dgx/miniconda3/lib/python3.7/site-packages/ase
    pandas     0.25.2   /home/dgx/miniconda3/lib/python3.7/site-packages/pandas
```

## How to use

The executable command of maptool is `mpt`,

```
mpt
```

you may select the corresponding function according to linux menu

```
+--------------------------------------------------------------------+
|                                     _              _               |
|               _ __ ___   __ _ _ __ | |_ ___   ___ | |              |
|              | '_ ` _ \ / _` | '_ \| __/ _ \ / _ \| |              |
|              | | | | | | (_| | |_) | || (_) | (_) | |              |
|              |_| |_| |_|\__,_| .__/ \__\___/ \___/|_|              |
|                              |_|                                   |
|                                                                    |
|                    0.1.dev52+g968793d.d20200420                    |
|                                                                    |
|                       Written by Wang haidi                        |
|                 URL  https://github.com/haidi-ustc                 |
|                Bug reports:(haidi@mail.ustc.edu.cn)                |
+--------------------------------------------------------------------+
======================== structural operation ========================
a1 >>> random operation
a2 >>> convert operation
a3 >>> build operation
a4 >>> cleave opeartion
a5 >>> strain operation
a6 >>> 2D structure operation
======================== structural analysis =========================
b1 >>> structure symmetry
b2 >>> structure finger print
b3 >>> structure difference
b4 >>> get primitive cell
b5 >>> get conventional cell
b6 >>> get XRD pattern
========================= vasp in/out tools ==========================
c1 >>> prepare input files
c2 >>> analysis output files
c3 >>> summary output files
====================== vasp calclation workflow=======================
d1 >>> optimize structure
d2 >>> calculate band structure
d3 >>> calculate band structure HSE06
d4 >>> calculate dos
d5 >>> calculate dos by HSE06
d6 >>> calculate elastic properties
d7 >>> calculate phonon
d8 >>> execute MD simulation
===================== Materials Project database =====================
e1 >>> get band/dos by mp-ID
e2 >>> get structure from materialsproject database
e3 >>> get properties by mp-ID
e4 >>> get phase graph
=========================== local database ===========================
f1 >>> check local database
f2 >>> get entry by l-ID
f3 >>> get entry by formula
f4 >>> get entry by element
f5 >>> insert entry into database
======================================================================


0  >>> back
88 >>> exit
your choce ?
--------------->>
```
For more information, please visit the 

## License
The code is open-source (licensed with a GPL license, see LICENSE)



