#!/usr/bin/env python
# coding: utf-8
import os
import itertools
import pandas as pd
from operator import itemgetter
from pymatgen.analysis.substrate_analyzer import SubstrateAnalyzer

try:
   from tqdm import tqdm
   TQDM=True
except ImportError:
   TQDM=False


def groupby_itemkey(iterable, item):
    """groupby keyed on (and pre-sorted by) itemgetter(item)."""
    itemkey = itemgetter(item)
    return itertools.groupby(sorted(iterable, key=itemkey), itemkey)

def match_substrate(film,substrates):
    all_matches = []
    sa = SubstrateAnalyzer()
    if TQDM:
        _substrates=tqdm(substrates)
    else:
        _substrates=substrate
    for s in _substrates:
        substrate = s["structure"]
        # Calculate all matches and group by substrate orientation
        matches_by_orient = groupby_itemkey(
            sa.calculate(film, substrate, lowest=True),
            "sub_miller")
        # Find the lowest area match for each substrate orientation
        lowest_matches = [min(g, key=itemgetter("match_area"))
                          for k, g in matches_by_orient]
        for match in lowest_matches:
            db_entry = {
                "sub_id": s["material_id"],
                "orient": " ".join(map(str, match["sub_miller"])),
                "sub_form": substrate.composition.reduced_formula,
                "film_orient": " ".join(map(str, match["film_miller"])),
                "area": match["match_area"],
            }
            all_matches.append(db_entry)
    df = pd.DataFrame(all_matches)
    df.set_index("sub_id", inplace=True)
    df.sort_values("area")
    return df
   
