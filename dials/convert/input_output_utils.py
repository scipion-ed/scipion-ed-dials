""" Functions for writing required json files (*.expt) for use with DIALS. """

import json
import msgpack


def writeJson(inputImages, fn='model.expt', idname="ExperimentList"):
    imageList = [img.clone() for img in inputImages]
    firstimage = imageList[0]
    lastimage = imageList[-1]
    templatepath = "{}/#####{}".format(firstimage.getDirName(),firstimage.getExtension())
    origin = [-firstimage.getBeamCenterMm()[0], firstimage.getBeamCenterMm()[1],
              -firstimage.getDistance()]
    exposure_time = []
    epoch = []
    for i in imageList:
        exposure_time.append(i.getExposureTime())
        epoch.append(i.getTwoTheta())
    beam = [
        {
            "direction": [
                0.0,
                0.0,
                1.0
            ],
            "wavelength": firstimage.getWavelength(),
            "divergence": 0.0,
            "sigma_divergence": 0.0,
            "polarization_normal": [
                0.0,
                1.0,
                0.0
            ],
            "polarization_fraction": 0.5,
            "flux": 0.0,
            "transmission": 1.0
        }
    ]
    detector = [
        {
            "panels": [
                {
                    "name": "Panel",
                    "type": "SENSOR_PAD",
                    "fast_axis": [
                        1.0,
                        0.0,
                        0.0
                    ],
                    "slow_axis": [
                        0.0,
                        -1.0,
                        0.0
                    ],
                    "origin": origin,
                    "raw_image_offset": [
                        0,
                        0
                    ],
                    "image_size":
                        firstimage.getDim(),
                    "pixel_size": [
                        firstimage.getPixelSize(),
                        firstimage.getPixelSize()
                    ],
                    "trusted_range": [
                        -1.0,
                        65535.0
                    ],
                    "thickness": 0.3,
                    "material": "Si",
                    "mu": 0.0,
                    "identifier": "",
                    "mask": [],
                    "gain": 1.0,
                    "pedestal": 0.0,
                    "px_mm_strategy": {
                        "type": "SimplePxMmStrategy"
                    }
                }
            ],
            "hierarchy": {
                "name": "",
                "type": "",
                "fast_axis": [
                    1.0,
                    0.0,
                    0.0
                ],
                "slow_axis": [
                    0.0,
                    1.0,
                    0.0
                ],
                "origin": [
                    0.0,
                    0.0,
                    0.0
                ],
                "raw_image_offset": [
                    0,
                    0
                ],
                "image_size": [
                    0,
                    0
                ],
                "pixel_size": [
                    0.0,
                    0.0
                ],
                "trusted_range": [
                    0.0,
                    0.0
                ],
                "thickness": 0.0,
                "material": "",
                "mu": 0.0,
                "identifier": "",
                "mask": [],
                "gain": 1.0,
                "pedestal": 0.0,
                "px_mm_strategy": {
                    "type": "SimplePxMmStrategy"
                },
                "children": [
                    {
                        "panel": 0
                    }
                ]
            }
        }
    ]

    goniometer = {
        "rotation_axis": firstimage.getRotationAxis(),
        "fixed_rotation": [
            1.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            1.0
        ],
        "setting_rotation": [
            1.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            1.0
        ],
    },
    scan = [
        {
            "image_range": [
                firstimage.getObjId(),
                lastimage.getObjId()
            ],
            "batch_offset": 0,
            "oscillation":
                firstimage.getOscillation(),
            "exposure_time": exposure_time,
            "epochs": epoch,
            "valid_image_ranges": {},
        }
    ]

    output = {
            "__id__": "{}".format(idname),
            "experiment": [
                {
                    "__id__": "Experiment",
                    "identifier": "",
                    "beam": 0,
                    "detector": 0,
                    "goniometer": 0,
                    "scan": 0,
                    "imageset": 0
                }
            ],
            "imageset": [
                {
                    "__id__": "ImageSequence",
                    "template": templatepath,
                    "mask": "",
                    "gain": "",
                    "pedestal": "",
                    "dx": "",
                    "dy": "",
                    "params": {
                        "dynamic_shadowing": "Auto",
                        "multi_panel": False
                },
            }
        ],
        "beam": beam,
        "detector": detector,
        "goniometer": goniometer,
        "scan": scan,
        "crystal": [],
        "profile": [],
        "scaling_model": []
    }

    with open(fn, 'w') as f:
        f.write(json.dumps(output, indent=4))
    return fn


def readRefl(reflFile, fn='reflections.txt', **kwargs):
    with open(reflFile, 'rb') as rf:
        contentsList = msgpack.unpackb(rf, use_bin_type=False, raw=True)
        for i in contentsList:
            print(type(i))
    # with open(fn, 'w') as f:
    #    f.write(contents)
    # print(contents)
    # return fn
