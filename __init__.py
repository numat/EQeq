"""
Python wrapper for EQeq, with a focus on streaming data.

This is intended 1. to allow users to automate + simplify their workflows and
2. to enable scaling of simulations to millions of structures.
"""
from ctypes import *

eqeq = cdll.LoadLibrary("./libeqeq.so")
eqeq.run.argtypes = (c_char_p, c_char_p)
eqeq.run.restype = c_char_p


def run(cif_structure, output_type="cif"):
    """Runs EQeq on the inputted structure, returning charge data.

    Args:
        cif_structure: Either a filename or data encoding a CIF file
        output_type: (Optional) Specifies the output type. Currently, options
            are "cif", "mol", "pdb", "car", and "files". The latter saves files
            of all possible output types.
    Returns:
        A string representing the charged crystal. Returns nothing if the
        output type is set to "files"
    """
    return eqeq.run(cif_structure, output_type.lower())


if __name__ == "__main__":
    from sys import argv
    print run(argv[1], argv[2])
