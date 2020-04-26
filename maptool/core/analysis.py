#!/usr/bin/env python3
# coding: utf-8

import numpy as np
from typing import List, Tuple
from nptyping import Array
from maptool import NAME
from maptool.util.utils import sepline,multi_structs
from maptool.io.read_structure import read_structures
from pymatgen import Structure,Molecule
from pymatgen.symmetry.analyzer import PointGroupAnalyzer,SpacegroupAnalyzer
from pymatgen.analysis.molecule_matcher import MoleculeMatcher
from pymatgen.analysis.diffraction.xrd import XRDCalculator


def structure_symmetry():
    structs,fnames=read_structures()
    multi_structs(structs,fnames)
    for struct,fname in zip(structs,fnames):
        if isinstance(struct,Structure):
            sa=SpacegroupAnalyzer(struct)
            print("{:<20} : {}".format('File Name',fname))
            print("{:<20} : {:<15}".format('Structure Type','periodicity'))
            print("{:<20} : {:<15}".format('Lattice Type',sa.get_lattice_type()))
            print("{:<20} : {:<15d}".format('Space Group ID',sa.get_space_group_number()))
            print("{:<20} : {:<15}".format('International Symbol',sa.get_space_group_symbol()))
            print("{:<20} : {:15}".format('Hall Symbol',sa.get_hall()))
            sepline()
        if isinstance(struct,Molecule):
            print("{:<20} : {}".format('File Name',fname))
            sa=PointGroupAnalyzer(struct)
            print("{:<20} : {:<15}".format('Structure Type','non-periodicity'))
            print("{:<20} : {:<15}".format('International Symbol',ast.get_pointgroup()))
    return True


def get_primitive_cell():
    structs,fnames=read_structures()
    multi_structs(structs,fnames)
    for struct,fname in zip(structs,fnames):
        sepline(ch='Primitive Cell',sp='-')
        ast=SpacegroupAnalyzer(struct)
        prim_st=ast.find_primitive()
        print(prim_st)
        sepline()
        print('save to '+NAME+'_primitive_'+fname+'.vasp')
        prim_st.to(filename=NAME+'_primitive_'+fname+'.vasp',fmt='poscar')
        sepline()
    return True


def get_conventional_cell():
    structs,fnames=read_structures()
    multi_structs(structs,fnames)
    for struct,fname in zip(structs,fnames):
        sepline(ch='conventional  cell',sp='-')
        ast=SpacegroupAnalyzer(struct)
        conv_st=ast.get_conventional_standard_structure()
        print(conv_st)
        sepline()
        print('save to '+NAME+'_convention_'+fname+'.vasp')
        conv_st.to(filename=NAME+'_conventional_'+fname+'.vasp',fmt='poscar')
        sepline()
    return True


def structure_finger_print():
    return None


def structures_difference(distance_tolerance=0.1,rcut=30):
    return None


def distance(struct1,struct2,rcut,pbc=False,):
    '''
    '''
    return None


def rdf(structures: List[Structure],
        r:                    float = 10,
        nbins:                  int = 80,
        range:  Tuple[float, float] = (0, 10),
        elem_pair:  Tuple[str, str] = ('', '')) -> Tuple[Array, Array]:
    '''
    Calculate radial distribution function for given trajectory

    @in
      - structures, [Structure], list of `Structure` or i.e. trajectory data
      - r, float, max raidus
      - nbins, int, number of bars in histogram
      - range, (float, float), plot range
      - elem_pair, (str, str), if the calculation of partial RDF is needed,
        pass this param. e.g. elem_pair = ('H', 'C'). Total RDF is returned if
        left empty.

    @out
      - (np.1darray, np.1darray), (rdf, radii) data.
    '''
    def shell_volumes(edges: Array):
        '''
        Calculates the shells' volume for given radii
        '''
        edges = edges ** 3
        return 4 * np.pi * (edges[1:] - edges[:-1]) / 3

    def rdf_helper(s:          Structure,
                   Aindices: Array[bool],
                   Bindices: Array[bool]) -> Array:
        '''
        The helper function. return the histogram data for each structure
        '''
        (centers, points, _, dist) = s.get_neighbor_list(r)
        validIndices = np.ones(dist.shape, dtype=bool)
        for (i, (cidx, pidx)) in enumerate(zip(centers, points)):
            validIndices[i] = Aindices[cidx] and Bindices[pidx]

        hist, _ = np.histogram(dist[validIndices], bins=nbins, range=range)
        return hist

    assert len(structures) > 0, 'Empty trajectory passed in'
    assert r > 0, 'Radius must be greater than 0'
    assert nbins > 0, 'nbin (number of bars of histogram) must be greater than 0'
    assert len(range) == 2 and\
        range[0] >= 0 and range[1] > 0 and\
        range[0] < range[1], 'Invalid range: range[0] and range[1] shoud in (0, ..] and range[0] < range[1]'
    st = structures[0]
    rho = st.num_sites / st.volume
    nsteps = len(structures)
    elem_array = np.array([s.species_string for s in st])
    if '' == elem_pair[0] and '' == elem_pair[1]:
        Aindices = np.ones(elem_array.shape, dtype=bool)
        Bindices = np.ones(elem_array.shape, dtype=bool)
    elif '' in elem_pair:
        raise Exception(f'Invalid element pair: {elem_pair}')
    elif elem_pair[0] not in elem_array or elem_pair[1] not in elem_array:
        raise Exception(f'Input elements f{elem_pair} not included in this structure')
    else:
        Aindices = elem_array == elem_pair[0]  # boolean array in order to index A
        Bindices = elem_array == elem_pair[1]  # boolean array in order to index B

    hist, edges = np.histogram([0], bins=nbins, range=range)
    hist *= 0

    for (i, s) in enumerate(structures):
        # print(f"Processing {i}")
        hist += rdf_helper(s, Aindices, Bindices)

    volumes = shell_volumes(edges)
    nA = np.sum(Aindices)  # number of atoms belong to element A
    nB = np.sum(Bindices)  # number of atoms belong to element B
    rdf = hist * st.num_sites / (volumes * nA * nB * rho * nsteps)
    radii = 0.5 * (edges[1:] + edges[:-1])
    return (rdf, radii)

def xrd(structure: Structure):
    '''
    Calculate XRD pattern of given structure, raw pattern data is returned

    @in
      - structure, Structure
    @out
      - x, np.1darray, theta
      - y, np.1darray, intensity
    '''
    c = XRDCalculator()
    pattern = c.get_pattern(structure)
    return (pattern.x, pattern.y)
