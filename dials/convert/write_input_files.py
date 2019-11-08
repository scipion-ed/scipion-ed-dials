import json


def writeJson(fn='model.expt', **kwargs):
    with open(fn, 'w') as f:
        f.write(json.dumps(kwargs, indent=4))
