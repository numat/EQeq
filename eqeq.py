"""
Python wrapper for EQeq, with a focus on streaming data.

This is intended 1. to allow users to automate + simplify their workflows and
2. to enable scaling of simulations to millions of structures.
"""
from ctypes import cdll, c_char_p, c_int, c_double, c_float
import os

ROOT = os.path.normpath(os.path.dirname(__file__))
eqeq = cdll.LoadLibrary(os.path.join(ROOT, "libeqeq.so"))
eqeq.run.argtypes = (c_char_p, c_char_p, c_double, c_float, c_int, c_char_p,
                     c_int, c_int, c_double, c_char_p, c_char_p)
eqeq.run.restype = c_char_p

DEFAULT_IONIZATION_PATH = os.path.join(ROOT, "ionizationdata.dat")
DEFAULT_CHARGE_PATH = os.path.join(ROOT, "chargecenters.dat")


def run(cif_structure, output_type="cif", l=1.2, h_i0=-2.0, charge_precision=3,
        method="direct", m_r=2, m_k=2, eta=50.0,
        ionization_data_path=DEFAULT_IONIZATION_PATH,
        charge_data_path=DEFAULT_CHARGE_PATH):
    """Runs EQeq on the inputted structure, returning charge data.

    Args:
        cif_structure: Either a filename or data encoding a CIF file.
        output_type: (Optional) Specifies the output type. Currently, options
            are "cif", "mol", "pdb", "car", and "files". The latter saves files
            of all possible output types.
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
    output_type, method = output_type.lower(), method.lower()
    if output_type not in ["cif", "pdb", "car", "mol", "json", "files"]:
        raise NotImplementedError("Output format '%s' is not supported!"
                                  % output_type)
    if method not in ["direct", "nonperiodic", "ewald"]:
        raise NotImplementedError("Method '%s' is not supported!" % method)

    return eqeq.run(cif_structure, output_type, l, h_i0, charge_precision,
                    method, m_r, m_k, eta, ionization_data_path,
                    charge_data_path)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Charge equilibration method "
                                     "for crystal structures.")
    parser.add_argument("input", type=str, help="An input cif file. Can be "
                        "either a filepath or the input data itself")
    parser.add_argument("--output-type", type=str, default="files",
                        help="Specifies the output type. Currently, options "
                        "are 'cif', 'mol', 'pdb', 'car', and 'files'. The "
                        "latter saves files of all possible output types")
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

    output = run(args.input, output_type=args.output_type, l=args.l,
                 h_i0=args.hi0, charge_precision=args.charge_precision,
                 method=args.method, m_r=args.mr, m_k=args.mk, eta=args.eta,
                 ionization_data_path=args.ionization_data_path,
                 charge_data_path=args.charge_data_path)
    print(output)
