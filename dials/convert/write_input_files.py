""" Functions for writing required json files (*.expt) for use with DIALS. """

import json


def writeJson(inputImages, fn='model.expt'):
    imageList = [img.clone() for img in inputImages]
    firstimage = imageList[0]
    lastimage = imageList[-1]
    templatepath = f"{firstimage.getDirName()}/#####{firstimage.getExtension()}"
    try:
        origin = [firstimage.getBeamCenter()[0], firstimage.getBeamCenter()[
            1], -firstimage.getDistance()]
    except Exception as e:
        print(e)
    exposure_time = []
    epoch = []
    for i in imageList:
        exposure_time.append(i.getExposureTime())
        epoch.append(i.getTwoTheta())
    beam = [{
        "direction": [
            0.0,
            0.0,
            1.0
        ],
        "transmission": 1.0,
        "polarization_normal": [
            0.0,
            1.0,
            0.0
        ],
        "divergence": 0.0,
        "polarization_fraction": 0.5,
        "flux": 0.0,
        "sigma_divergence": 0.0,
        "wavelength": firstimage.getWavelength()
    }]
    detector = [
        {
            "hierarchy": {
                "origin": [
                    0.0,
                    0.0,
                    0.0
                ],
                "fast_axis": [
                    1.0,
                    0.0,
                    0.0
                ],
                "name": "",
                "raw_image_offset": [
                    0,
                    0
                ],
                "slow_axis": [
                    0.0,
                    1.0,
                    0.0
                ],
                "material": "",
                "mask": [],
                "thickness": 0.0,
                "mu": 0.0,
                "gain": 1.0,
                "trusted_range": [
                    0.0,
                    0.0
                ],
                "image_size": [
                    0,
                    0
                ],
                "px_mm_strategy": {
                    "type": "SimplePxMmStrategy"
                },
                "pedestal": 0.0,
                "identifier": "",
                "type": "",
                "children": [
                    {
                        "panel": 0
                    }
                ],
                "pixel_size": [
                    0.0,
                    0.0
                ]
            },
            "panels": [
                {
                    "origin": origin,
                    "fast_axis": [
                        1.0,
                        0.0,
                        0.0
                    ],
                    "name": "Panel",
                    "raw_image_offset": [
                        0,
                        0
                    ],
                    "slow_axis": [
                        0.0,
                        -1.0,
                        0.0
                    ],
                    "material": "Si",
                    "mask": [],
                    "thickness": 0.3,
                    "mu": 0.0,
                    "gain": 1.0,
                    "trusted_range": [
                        -1.0,
                        65535.0
                    ],
                    "image_size": [
                        0,
                        0
                    ],
                    "px_mm_strategy": {
                        "type": "SimplePxMmStrategy"
                    },
                    "pedestal": 0.0,
                    "identifier": "",
                    "type": "SENSOR_PAD",
                    "pixel_size": [
                        firstimage.getPixelSize(),
                        firstimage.getPixelSize()
                    ]
                }
            ]
        }
    ]
    goniometer = [
        {
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
            "rotation_axis": [
                -0.6203929740893525,
                -0.7842911179533836,
                0.0
            ]
        }
    ],
    scan = [
        {
            "exposure_time": exposure_time,
            "batch_offset": 0,
            "oscillation": [
                firstimage.getOscillation()
            ],
            "valid_image_ranges": {},
            "epochs": epoch,
            "image_range": [
                firstimage.getIndex(),
                lastimage.getIndex()
            ],
        }
    ]

    output = [
        {
            "__id__": "DataBlock",
            "imageset": [{
                "__id__": "ImageSweep",
                "template": templatepath,
                "mask": None,
                "gain": None,
                "pedestal": None,
                "dx": None,
                "dy": None,
                "beam": 0,
                "detector": 0,
                "goniometer": 0,
                "scan": 0,
                "params": {
                    "dynamic_shadowing": "Auto",
                    "multi_panel": False
                },
            }],
            "beam": beam,
            "detector": detector,
            "goniometer": goniometer,
            "scan": scan,
        }
    ]

    with open(fn, 'w') as f:
        f.write(json.dumps(output, indent=4))
    return fn
