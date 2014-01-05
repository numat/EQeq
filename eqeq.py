"""
Python wrapper for EQeq, with a focus on streaming data.

This is intended 1. to allow users to automate + simplify their workflows and
2. to enable scaling of simulations to millions of structures.
"""
from ctypes import cdll, c_char_p, c_int, c_double, c_float
import json_formatter as json
import os
import sys

try:
    import format_converter
except:
    sys.stderr.write("Openbabel not found. Format conversion unsupported.\n"
                     "Only use cif files for output, and cif/mol/pdb/car/files"
                     " for output.\n")

sys.stderr.write("This version of EQeq is deprecated. Please switch to the "
                 "version maintained by the openbabel project (see the README "
                 "for more).\n")

eqeq = cdll.LoadLibrary("/usr/lib/libeqeq.so")
eqeq.run.argtypes = (c_char_p, c_char_p, c_double, c_float, c_int, c_char_p,
                     c_int, c_int, c_double, c_char_p, c_char_p)
eqeq.run.restype = c_char_p

ROOT = os.path.normpath(os.path.dirname(__file__))
DEFAULT_IONIZATION_PATH = os.path.join(ROOT, "ionizationdata.dat")
DEFAULT_CHARGE_PATH = os.path.join(ROOT, "chargecenters.dat")


def run(structure, input_type="cif", output_type="cif", l=1.2, h_i0=-2.0,
        charge_precision=3, method="ewald", m_r=2, m_k=2, eta=50.0,
        ionization_data_path=DEFAULT_IONIZATION_PATH,
        charge_data_path=DEFAULT_CHARGE_PATH):
    """Runs EQeq on the inputted structure, returning charge data.

    Args:
        structure: Either a filename or data encoding a chemical.
        input_type: (Optional) Specifies input type. Can be anything supported
            by openbabel, as well as "json"
        output_type: (Optional) Specifies the output type. Currently, options
            are "cif", "mol", "pdb", "car", "json", "list", and "files". The
            first four return modified chemical data formats, "list" returns
            a Python object, "json" is that object serialized, and "files"
            saves files of all possible output types.
        l: (Optional) Lambda, the dielectric screening parameter.
        h_i0: (Optional) The electron affinity of hydrogen.
        charge_precision: (Optional) Number of decimals to use for charges.
        method: (Optional) Method to use. Can be "direct" (default),
            "nonperiodic", or "ewald".
        m_r: (Optional) Number of unit cells to consider in "real space". This
            is measured radially, so m_r = 1 evaluates 27 unit cells.
        m_k: (Optional) Number of unit cells to consider in "frequency space".
            This is measured radially, so m_k = 1 evaluates 27 unit cells.
        eta: (Optional) Ewald splitting parameter
        ionization_data_path: (Optional) A path to the file containing ion-
            ization data. By default, assumes the data is in the EQeq folder
            and saved as "ionizationdata.dat".
        charge_data_path: (Optional) A path to the file containing charge-
            center data. By default, assumes the data is in the EQeq folder
            and saved as "chargecenters.dat".
    Returns:
        A string representing the charged crystal. Returns nothing if the
        output type is set to "files"
    """
    # Error handling on string params. Should spare users some annoyance.
    o, m = output_type.lower(), method.lower()
    if o not in ["cif", "pdb", "car", "mol", "json", "list", "files"]:
        raise NotImplementedError("Output format '%s' is not supported!" % o)
    if m not in ["direct", "nonperiodic", "ewald"]:
        raise NotImplementedError("Method '%s' is not supported!" % m)
    # If linked to openbabel, use it to handle json interconversion externally
    if input_type != "cif":
        structure = format_converter.convert(structure, input_type, "cif")
    structure = structure.replace("\t", "  ")
    # Calls libeqeq.so's run method, returning a string of data
    result = eqeq.run(structure, ("json" if output_type == "list" else
                      output_type), l, h_i0, charge_precision, method, m_r,
                      m_k, eta, ionization_data_path, charge_data_path)
    if output_type == "list":
        return json.loads(result)
    # This option appends atoms in json/object data with a "charge" attribute
    if output_type == "json":
        obj = format_converter.convert(structure, "cif", "object")
        result = json.loads(result)
        for atom, charge in zip(obj["atoms"], result):
            atom["charge"] = charge
        result = json.dumps(obj)
    return result


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Charge equilibration method "
                                     "for crystal structures.")
    parser.add_argument("input", type=str, help="An input cif file. Can be "
                        "either a filepath or the input data itself")
    parser.add_argument("--input-type", type=str, default="cif",
                        help="Specifies input type. Can be anything supported "
                        "by openbabel and 'json'")
    parser.add_argument("--output-type", type=str, default="files",
                        help="Specifies the output type. Currently, options "
                        "are 'cif', 'mol', 'pdb', 'car', 'json', and 'files'. "
                        "The latter saves files of all possible output types")
    parser.add_argument("--l", type=float, default=1.2,
                        help="Lambda, the dielectric screening parameter")
    parser.add_argument("--hi0", type=float, default=-2.0,
                        help="The electron affinity of hydrogen")
    parser.add_argument("--charge-precision", type=int, default=3,
                        help="Number of decimals to use for charges")
    parser.add_argument("--method", type=str, default="direct",
                        help="Method to use. Can be 'direct' (default), "
                        "'nonperiodic', or 'ewald'")
    parser.add_argument("--mr", type=int, default=2, help="Number of unit "
                        "cells to consider in 'real space'. This is measured "
                        "radially, so m_r=1 evaluates 27 unit cells")
    parser.add_argument("--mk", type=int, default=2, help="Number of unit "
                        "cells to consider in 'frequency space'. This is "
                        "measured radially, so m_k=1 evaluates 27 unit cells")
    parser.add_argument("--eta", type=float, default=50.0,
                        help="Ewald splitting parameter")
    parser.add_argument("--ionization-data-path", type=str,
                        default=DEFAULT_IONIZATION_PATH, help="A path to the "
                        "file containing ionization data. By default, assumes "
                        "the data is in the EQeq folder and saved as "
                        "'ionizationdata.dat'")
    parser.add_argument("--charge-data-path", type=str,
                        default=DEFAULT_CHARGE_PATH, help="A path to the file "
                        "containing charge-center data. By default, assumes "
                        "the data is in the EQeq folder and saved as "
                        "'chargecenters.dat'")
    args = parser.parse_args()

    output = run(args.input, input_type=args.input_type,
                 output_type=args.output_type, l=args.l,
                 h_i0=args.hi0, charge_precision=args.charge_precision,
                 method=args.method, m_r=args.mr, m_k=args.mk, eta=args.eta,
                 ionization_data_path=args.ionization_data_path,
                 charge_data_path=args.charge_data_path)
    print(output)
