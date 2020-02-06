""" Functions for writing required json files (*.expt) for use with DIALS. """

import json
import msgpack
import numpy as np


def writeJson(inputImages, fn='model.expt', idname="ExperimentList"):
    imageList = [img.clone() for img in inputImages]
    firstimage = imageList[0]
    lastimage = imageList[-1]
    templatepath = "{}/#####{}".format(firstimage.getDirName(),
                                       firstimage.getExtension())
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


def readRefl(reflFile, **kwargs):
    with open(reflFile, 'rb') as f:
        buf = msgpack.unpack(f)

    reflFileIdentifier = buf[0].decode()
    version = buf[1]
    nrows = buf[2][b'nrows']
    identifier_dict = buf[2][b'identifiers']
    data_dict = buf[2][b'data']
    data = {}
    for k, v in data_dict.items():
        data[k.decode()] = np.array(extractRefls(v))
    return reflFileIdentifier, version, nrows, identifier_dict, data


def extractRefls(v):
    dtype = v[0]
    size = v[1][0]
    if dtype == b"int":
        data = np.frombuffer(v[1][1], dtype=np.int32)
        assert len(data) == size
    elif dtype == b"int6":
        data = np.frombuffer(v[1][1], dtype=np.int32)
        assert len(data) == size * 6
        data = data.reshape((size, 6))
    elif dtype == b"std::size_t":
        data = np.frombuffer(v[1][1], dtype=np.uint64)
        assert len(data) == size
    elif dtype == b"double":
        data = np.frombuffer(v[1][1], dtype=np.float64)
        assert len(data) == size
    elif dtype == b"vec3<double>":
        data = np.frombuffer(v[1][1], dtype=np.float64)
        assert len(data) == size * 3
        data = data.reshape((size, 3))
    else:
        data = None
    return data


def writeRefl(inputRefls, fn='reflections.refl', reflFileIdentifier='dials::af::reflection_table', version=1, nrows=None, identifier_dict={}, **kwargs):
    # TODO: Get all data from database and add additional ones
    if nrows is not None:
        nrows = nrows
    else:
        nrows = inputRefls.getSpots()
    
    # TODO: Do the inverse conversion of extractRefls above
    # TODO: Create the overall structure to pack
    # Initialise the parts of the structure to pack
    data = {'bbox': ['int6', [nrows, None]],
            'flags': ['std::size_t', [nrows, None]],
            'id': ['int', [nrows, None]],
            'intensity.sum.value': ['double', [nrows, None]],
            'intensity.sum.variance': ['double', [nrows, None]],
            'n_signal': ['int', [nrows, None]],
            'panel': ['std::size_t', [nrows, None]],
            'shoebox': ['Shoebox<>', [nrows, None]],
            'xyzobs.px.value': ['vec3<double>', [nrows, None]],
            'xyzobs.px.variance': ['vec3<double>', [nrows, None]]
            }
    content_dict = {'identifiers': identifier_dict,
                    'nrows': nrows, 'data': data}

    # TODO: Pack all data using msgpack
    output = [reflFileIdentifier, version, content_dict]
    with open(fn, 'wb') as f:
        f.write(msgpack.packb(output))
    return fn
