## Introduction
Maptool (Material Project Tool) is an open source python package for pre- and post- process input and output of first principles simulation software, which is based on Pymatgen and ASE and other open source python packages. Currently, the VASP code are mainly supported, later it will support more DFT codes. The package allows for generating, building, modifying and analysing of both molecular and crystal structure. Besides, maptool can automatically generating the input files for different DFT code and analysing the output, all results are saved with pretty format and figure based on matplotlib. See the [documentation](https://github.com/haidi-ustc/maptool/tree/master/doc/) for information about installation and usage.

# Current Features:
* Random generation of atomic (3D, 2D, and 1D) crystals and point group clusters for a given symmetry group and stoichiometry
* Building of cluster, shell strcuture and surface
* Add strain to crystal structure
* Analysing of symmetry of molecule and crystal structure
* Structure transformation to cif, xsf, xyz, mol or POSCAR files via pymatgen and ASE
* General operation for 2D materials (splitting, strain, bending and so on)
* VASP input generation and output analysing and visualization
* Online retrieval function based on pymatgen
* Local database for storing the previous result

## Dependencies:
* [SciPy 1.3.1](https://www.scipy.org/install.html)
* [NumPy 1.17.3](https://www.scipy.org/scipylib/download.html)
* [Pandas 0.25.2](https://pandas.pydata.org/getpandas.html)
* [Pymatgen 2019.10.16](http://pymatgen.org/#getting-pymatgen)
* [ASE 3.19.0](https://pypi.org/project/ase)

### Optional:
* [openbabel 2.4.1 (Python bindings)](http://openbabel.org/wiki/Main_Page) (allows for additional molecule file formats)
* [pymongo](https://api.mongodb.com/python) (allows for database operation)

## Installation

```
$ git clone git@github.com:haidi-ustc/maptool.git
$ cd maptool
$ python setup.py install
```
