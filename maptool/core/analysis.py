#!/usr/bin/env python3
# coding: utf-8

import numpy as np
from typing import List, Tuple
from nptyping import Array
from maptool import NAME
from maptool.util.utils import (
    sepline,
    multi_structs
)
from maptool.io.read_structure import read_structures
from pymatgen import (
    Structure,
    Molecule
)
from pymatgen.symmetry.analyzer import (
    PointGroupAnalyzer,
    SpacegroupAnalyzer
)
from pymatgen.analysis.molecule_matcher import MoleculeMatcher
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.analysis.structure_prediction.volume_predictor import DLSVolumePredictor
import matplotlib.pyplot as plt


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
            pga=PointGroupAnalyzer(struct)
            print("{:<20} : {:<15}".format('Structure Type','non-periodicity'))
            print("{:<20} : {:<15}".format('International Symbol',str(pga.get_pointgroup())))
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

def xrd(structure:       Structure,
        two_theta_range: Tuple[int, int] = (0, 120),
        sigma:                     float = 0.05,
        nsample:                     int = 5000,
        fig_name:                    str = "XRD.png",
        peak_raw_fname:              str = "",
        plot_dat_fname:              str = ""):
    '''
    Calculate XRD pattern of given structure, raw pattern data is returned

    @in
      - structure, Structure
      - two_theta_range, (int, int), range of 2Theta in XRD pattern
      - sigma, float, gaussian smearing width
      - nsample, int, how many points in final plot
      - fig_name, str, name of the figure file to be written,
        default is XRD.png
      - peak_raw_fname, str, name of the raw data file containing 2Theta and
        intensity.
      - plot_dat_fname, str, name of the plot data file. That plot data was
        used to plot the XRD pattern.
    @out
      - x, np.1darray, theta
      - y, np.1darray, intensity
      - plot_x, np.1darray, theta of final data in plotting
      - plot_y, np.1darray, intensity of final data in plotting
    '''
    def _smearing(x, y,
                  sigma=sigma,
                  nsample=nsample) -> Tuple[Array, Array]:
        '''
        Apply smearing to discret delta functions to get a continuous array
        @in
          - x: x coordinates of peak
          - y: peak heights
        @out
          - res_x
          - res_y
        '''
        def _smearing_helper(x, x0, sigma=sigma):
            return np.exp(-(x - x0)**2 / (2 * sigma**2))
        assert x.shape == y.shape
        npeaks = x.shape[0]
        res_x = np.linspace(two_theta_range[0],
                            two_theta_range[1], num=nsample)
        smear_res = np.empty((npeaks, nsample))
        for ipeak in range(npeaks):
            smear_res[ipeak] = _smearing_helper(res_x, x[ipeak]) * y[ipeak]
            pass
        res_y = np.sum(smear_res, axis=(0))
        return (res_x, res_y)

    c = XRDCalculator()
    p = c.get_pattern(structure,
                      two_theta_range=two_theta_range)
    (x, y) = _smearing(p.x, p.y, sigma=sigma)
    plt.figure()
    plt.plot(x, y)
    plt.vlines(p.x, ymin=-5, ymax=-1)
    plt.ylim(-5, 110)
    plt.xlabel(r"$2\theta$ ($^\circ$)")
    plt.ylabel(r"Intensities (scaled)")
    if "" != fig_name:
        plt.savefig(fig_name, dpi=800, linewidth=0.01)

    if "" != peak_raw_fname:
        with open(peak_raw_fname, 'w') as f:
            to_be_written = "# {:^10} {:^12} {:^9}\n".format(
                "2Theta", "Intensity", "Miller")
            for (_x, _y, hkls) in zip(p.x, p.y, p.hkls):
                label = ",".join([str(hkl['hkl']) for hkl in hkls])
                to_be_written += " {:11.7f} {:11.7f} {}\n".format(_x, _y, label)
            print(to_be_written, file=f)

    if "" != plot_dat_fname:
        with open(plot_dat_fname, 'w') as f:
            to_be_written = "# {:^10} {:^12}\n".format(
                "2Theta", "Intensity")
            for (_x, _y) in np.vstack([x, y]).T:
                to_be_written += " {:11.7f} {:11.7f}\n".format(_x, _y)
            print(to_be_written, file=f)
    return (p.x, p.y, x, y)


def structure_dedup(structures: List[Structure],
                    fnames:     List[str] = []) -> Tuple[List[Structure],
                                                         List[str]]:
    '''
    Deduplicate the given structures via pymatgen.analysis.structure_matcher
    @in
      - structures, [Structure], given structures to be deduplicate
      - fnames, [str], file names corresponding to structures. If left empty,
        an empty list will be returned
    @out
      - [Structure], deduplicated structures

    '''
    sm = StructureMatcher()

    def _match(st_ref: Structure,
               st_lst: List[Structure]):
        for st in st_lst:
            if sm.fit(st, st_ref, symmetric=True):
                return True
        return False

    assert len(fnames) == 0 or len(fnames) == len(structures),\
        'fname list should be empty or have same size with structures list'

    if len(structures) == 0 or len(structures) == 1:
        return structures, fnames
    else:
        clist = []
        ilist = []
        for (i, st) in enumerate(structures):
            if not _match(st, clist):
                clist.append(st.copy())
                ilist.append(i)
        flist = []
        if len(fnames) != 0:
            flist = [fnames[i] for i in ilist]
        return clist, flist


def volume_predict(structure: Structure) -> (
        Structure, float):
    '''
    Returns the scaled structure and volume of given structure according to
    pymatgen's prediction

    @in
      - structure, Structure, given structure
    @out
      - Structure, scaled structure
      - float, volume of the scaled strcture
    '''
    dls = DLSVolumePredictor()
    st = structure
    dls_st = dls.get_predicted_structure(st)
    volume = dls.predict(st)
    return (dls_st, volume)
