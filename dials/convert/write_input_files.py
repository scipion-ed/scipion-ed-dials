""" Functions for writing required json files (*.expt) for use with DIALS.
Might be more suitable to perform with dials or dxtbx as a module. """

import json


def writeJson(fn='model.expt', **kwargs):
    with open(fn, 'w') as f:
        f.write(json.dumps(kwargs, indent=4))
