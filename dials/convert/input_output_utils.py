""" Functions for writing required json files (*.expt) for use with DIALS. """

import os
import json
import msgpack
import numpy as np
import shutil


def writeJson(inputImages, fn='model.expt', idname="ExperimentList", overwriteModel=False):
    if overwriteModel is False and os.path.exists(fn):
        return
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


def readRefl(reflFile, fn='reflections.txt', **kwargs):
    with open(reflFile, 'rb') as f:
        buf = msgpack.unpack(f, strict_map_key=False)

    reflFileIdentifier = buf[0]
    version = buf[1]
    nrows = buf[2]['nrows']
    identifier_dict = buf[2]['identifiers']
    data_dict = buf[2]['data']
    data = {}
    for k, v in data_dict.items():
        data[k] = np.array(extractRefls(v))
    return reflFileIdentifier, version, nrows, identifier_dict, data


def extractRefls(v):
    dtype = v[0]
    size = v[1][0]
    if dtype == "int":
        data = np.frombuffer(v[1][1], dtype=np.int32)
        assert len(data) == size
    elif dtype == "int6":
        data = np.frombuffer(v[1][1], dtype=np.int32)
        assert len(data) == size * 6
        data = data.reshape((size, 6))
    elif dtype == "std::size_t":
        data = np.frombuffer(v[1][1], dtype=np.uint64)
        assert len(data) == size
    elif dtype == "double":
        data = np.frombuffer(v[1][1], dtype=np.float64)
        assert len(data) == size
    elif dtype == "vec3<double>":
        data = np.frombuffer(v[1][1], dtype=np.float64)
        assert len(data) == size * 3
        data = data.reshape((size, 3))
    else:
        data = None
    return data


def compressRefl(data):
    data_dict = {}
    return data_dict


def writeRefl(inputSpots, fn='reflections.refl', **kwargs):
    # Workaround to use existing file
    if type(inputSpots) is str:
        try:
            shutil.copy(inputSpots, fn)
        except shutil.SameFileError:
            pass

            # FIXME: Make the below implementation work
    else:
        # FIXME: Convert to use SetOfSpots
        [reflFileIdentifier, version, nrows, identifier_dict, data] = inputSpots

        data_dict = compressRefl(data)

        header_dict = {
            'nrows': nrows,
            'identifiers': identifier_dict,
            'data': data_dict
        }
        output = [reflFileIdentifier.encode(), version, header_dict]

        with open(fn, 'wb') as f:
            f.write(msgpack.packb(output))


def writeRefinementPhil(fn='refinement.phil', **kwargs):
    template = ["refinement {",
                "   parameterisation {",
                "        beam {",
                "            fix = all *in_spindle_plane out_spindle_plane *wavelength",
                "            }",
                "        crystal {",
                "            fix = all cell orientation",
                "            }",
                "        detector {",
                "            fix = all position orientation distance",
                "            }",
                "        goniometer {",
                "            fix = *all in_beam_plane out_beam_plane",
                "      }",
                "   }",
                "    reflections {",
                "    outlier {",
                "      algorithm = null *auto mcd tukey sauter_poon",
                "    }",
                "  }",
                "}"]

    with open(fn, 'w') as f:
        f.write("\n".join(template))


def copyDialsFile(originalDialsFile, fn=None):
    try:
        shutil.copy(originalDialsFile, fn)
    except shutil.SameFileError:
        pass
