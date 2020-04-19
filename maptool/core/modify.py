#!/usr/bin/env python3

import os
import random
from random import randint
import itertools
import numpy as np
from typing import List, Tuple, Dict
from pymatgen import (
  Structure,
  Molecule,
  Lattice,
)
from pymatgen.core.periodic_table import is_valid_symbol
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.surface import (
  generate_all_slabs,
  SlabGenerator,
  get_symmetrically_distinct_miller_indices,
  )
from pymatgen.transformations.site_transformations \
  import ReplaceSiteSpeciesTransformation

class StructureChanger:
  def __init__(self,
               structure: Structure):
    self.old_structure = structure
    self.operations = []

  def scale_volume(self,
                   scale: float = 0.01) -> Structure:
    '''
    Scale the lattice with given parameter.
    Vol_new = Vol_old * (1 + scale)

    @in
      - scale, float

    @out
      - Structure
    '''
    struct = self.old_structure.copy()
    struct.scale_lattice(
      (1 + scale) * self.old_structure.volume
    )
    self.operations.append({'V scaling': scale})
    return struct

  def apply_defect(self,
                   number:          int,
                   defect_type:     str = 'vac',
                   supercell: List[int] = [1, 1, 1],
                   element:         str = '') -> Structure:
    '''
    Apply the defects, i.e. make changes to the structure including
    removing, substitution and etc.

    @in
      - defect_type, str:
        - vac: randomly remove some sites with given number of atoms
        - subs: randomly substitute some sites with given element
                a number of substitution operations is need.
        - inte: not implemented yet
        - vint: not implemented yet
      - supercell, list of int, len == 3

    @out
      - Structure
    '''
    assert number >= 0, f"Number of defects should be non-negative"
    assert len(supercell) == 3, f"Length of supercell list must be 3"

    sc = self.old_structure.copy()
    sc.make_supercell(supercell)
    indices = random.sample(range(len(sc)), number)  # Select the apply site
    if 'vac' == defect_type:
      rest_id = [i for i in range(len(sc)) if i not in indices]
      new_structure = sc[rest_id].copy()
      self.operations.append({'vacancy': number})
      return new_structure

    elif 'subs' == defect_type:
      assert is_valid_symbol(element), f'Invalid element input: "{element}"'
      ops = {i: element for i in indices}
      trans = ReplaceSiteSpeciesTransformation(ops)
      new_structure = trans.apply_transformation(sc)
      new_structure.sort()
      self.operations.append({'substitution': number})
      return new_structure

    elif 'inte' == defect_type:
      raise Exception(f'This type of defect ({defect_type}) is still in developing')

    elif 'vint' == defect_type:
      raise Exception(f'This type of defect ({defect_type}) is still in developing')

    else:
      raise Exception(f'Invalid defect type input: "{defect_type}"')

  def get_slabs(self,
                max_index:                      int = 1,
                min_slab_size:                float = 5.0,
                min_vacuum_size:              float = 15.0,
                fix_vacuum_size:               bool = False,
                bonds: Dict[Tuple[str, str], float] = None,
                tolerance:                    float = 0.001,
                max_broken_bonds:               int = 0,
                is_lll_reduce:                 bool = False,
                is_center_slab:                bool = False,
                is_primitive:                  bool = True,
                max_normal_search:              int = None,
                is_symmetrize:                 bool = False,
                is_repair:                     bool = False,
                is_in_unit_planes:             bool = False) -> List[Structure]:
    '''
    Search and return the slab found in the given structure.

    '''

    st = self.old_structure.copy()
    all_slabs = []
    for miller in get_symmetrically_distinct_miller_indices(st, max_index):
        if fix_vacuum_size:
            pass

        else:
            vacuum_size = random.random() * min_vacuum_size

        gen = SlabGenerator(st,
                            miller,
                            min_slab_size,
                            vacuum_size,
                            lll_reduce=is_lll_reduce,
                            center_slab=is_center_slab,
                            primitive=is_primitive,
                            max_normal_search=max_normal_search,
                            in_unit_planes=is_in_unit_planes)
        slabs = gen.get_slabs(bonds=bonds,
                              tol=tolerance,
                              symmetrize=is_symmetrize,
                              max_broken_bonds=max_broken_bonds,
                              repair=is_repair)

        if len(slabs) > 0:
            all_slabs.extend(slabs)
#                for slab in all_slabs:
#                    slab.make_supercell(supercell)

    self.operations.append({'slabs': len(slabs)})
    return all_slabs
