""" Functions related to writing and handling phil files """


import re


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


def writeRestraintsPhil(fn='restraints.phil', values=None, sigmas=None):

    # Helper function to sanitise the input
    def fixInput(inputString):
        if inputString is None:
            return None
        # Change from string splitting numbers with any combination of
        # commas and white spaces into a list of floats
        inputList = [float(s) for s in re.split(
            '[,\s]', inputString.strip()) if s]

        # Ensure that the list has 6 elements by adding or removing the last entries
        while len(inputList) < 6:
            inputList.extend([inputList[-1]])
        del inputList[6:]

        # Convert the list of floats back into a string with comma as separator
        fixedString = ','.join([str(s) for s in inputList])
        return fixedString

    # No point in writing the file if there are no values to tie to
    if values is None:
        return

    contentString = [
        f"                        values={fixInput(values)}"
    ]
    # Only add the string with sigmas if there are any specific sigmas to add
    if sigmas is not None:
        contentString += [
            f"                        sigmas={fixInput(sigmas)}"]
    finString = '\n'.join(contentString)
    template = [
        "refinement",
        "{",
        "    parameterisation",
        "    {",
        "        crystal",
        "        {",
        "            unit_cell",
        "            {",
        "                restraints",
        "                {",
        "                    tie_to_target",
        "                    {",
        f"{finString}",
        "                    }",
        "                }",
        "            }",
        "        }",
        "    }",
        "}",
    ]

    with open(fn, 'w') as f:
        f.write("\n".join(template))
    return fn
