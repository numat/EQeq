"""
Methods to interconvert between json and other (cif, mol, smi, etc.) files
"""

import pybel
ob = pybel.ob
import numpy as np
from numpy.linalg import norm, inv

import json_formatter as json

table = ob.OBElementTable()


def convert(data, in_format, out_format, pretty=False):
    """Converts between two inputted chemical formats.

    Args:
        data: A string representing the chemical file to be converted. If the
            `in_format` is "json", this can also be a Python object
        in_format: The format of the `data` string. Can be "json" or any format
            recognized by Open Babel
        out_format: The format to convert to. Can be "json" or any format
            recognized by Open Babel
        pretty: (Optional) If True and `out_format` is "json", will pretty-
            print the output for human readability
    Returns:
        A string representing the inputted `data` in the specified `out_format`
    """

    # Decide on a json formatter depending on desired prettiness
    dumps = json.dumps if pretty else json.compress

    # If it's a json string, load it
    if in_format == "json" and isinstance(data, basestring):
        data = json.loads(data)

    # A little "hack" to format inputted json
    if in_format == "json" and out_format == "json":
        return json.dumps(data)

    # These are converted manually to retain crystallographic information
    if in_format == "json" and out_format == "cif":
        return json_to_cif(data)

    # These use the open babel library to interconvert, with additions for json
    mol = (json_to_pybel(data) if in_format == "json" else
           pybel.readstring(in_format.encode("ascii"),
                            "".join(i for i in data if ord(i) < 128)
                            .encode("ascii")))

    # Infer structure in cases where the input format has no specification
    if not mol.OBMol.HasNonZeroCoords():
        mol.make3D()
    mol.OBMol.Center()

    # EQeq takes a specific cif format that openbabel does not output.
    # This manually overrides that.
    if out_format == "cif":
        return json_to_cif(pybel_to_json(mol))

    if out_format == "object":
        return pybel_to_json(mol)
    elif out_format == "json":
        return dumps(pybel_to_json(mol))
    else:
        return mol.write(out_format)


def json_to_cif(mof):
    """Converts python mof data structure to a cif file.

    Args:
        mof: A Python data structure analogous to `structures.Mof`
    Returns:
        A cif file containing atom and periodic connection data
    """
    va, vb, vc = [np.array(v) for v in mof["periodic_connections"]]
    a, b, c = [norm(i) for i in [va, vb, vc]]
    alpha, beta, gamma = [get_angle(i, j) * 180 / np.pi
                          for i, j in [[vb, vc], [va, vc], [vb, va]]]

    # Headers and unit cell information (note the space groups are wrong)\n"
    cif = ("data_functionalizedCrystal\n"
           " _audit_creation_method\t'MofGen v2!'\n"
           "_symmetry_space_group_name_H-M\t'P1'\n"
           "_symmetry_Int_Tables_number\t1\n"
           "_symmetry_cell_setting\ttriclinic\n"
           "loop_\n"
           "_symmetry_equiv_pos_as_xyz\n"
           "   x,y,z\n"
           "_cell_length_a\t%f\n"
           "_cell_length_b\t%f\n"
           "_cell_length_c\t%f\n"
           "_cell_angle_alpha\t%f\n"
           "_cell_angle_beta\t%f\n"
           "_cell_angle_gamma\t%f\n"
           "loop_\n"
           "_atom_site_label\n"
           "_atom_site_type_symbol\n"
           "_atom_site_fract_x\n"
           "_atom_site_fract_y\n"
           "_atom_site_fract_z\n") % (a, b, c, alpha, beta, gamma)

    # Adding atom positional data in fractional coordinates
    if "building_blocks" in mof:
        mof["atoms"] = [atom for bb in mof["building_blocks"]
                        for atom in bb["atoms"]]
    cif += "\n".join("%s%d\t%s\t%.5f\t%.5f\t%.5f"
                     % tuple([atom["element"], i, atom["element"]] +
                     cartesian_to_fractional(atom["location"], va, vb, vc))
                     for i, atom in enumerate(mof["atoms"]))
    return cif + "\n"


def get_angle(first_vector, second_vector):
    """Returns the angle, in radians, between two vectors."""
    return np.arccos(np.clip(np.dot(first_vector, second_vector) /
                    (norm(first_vector) * norm(second_vector)), -1, 1))


def cartesian_to_fractional(location, va, vb, vc):
    """Convert cartesian to fractional coordinates using unit cell vectors.

    Args:
        location: The cartesian coordinates to be converted
        va, vb, vc: The three crystal vectors of the structure
    Returns:
        A three-element list containing the fractional coordinates
    """
    # `cartesian` = matrix(`vectors`) * `fractional`. Invert to solve.
    frac = inv(np.matrix([va, vb, vc]).transpose()).dot(np.array(location))
    # This constrains the fractional coordinates to [0, 1]
    return [f % 1 for f in np.nditer(frac)]


def json_to_pybel(data, center=True):
    """Converts python data structure to pybel.Molecule.

    This will infer bond data if not specified.

    Args:
        data: The loaded json data of a molecule, as a Python object
        center: (Optional) Centers the coordinates of the outputted molecule
    Returns:
        An instance of `pybel.Molecule`
    """
    obmol = ob.OBMol()
    obmol.BeginModify()
    for atom in data["atoms"]:
        obatom = obmol.NewAtom()
        obatom.SetAtomicNum(table.GetAtomicNum(str(atom["element"])))
        obatom.SetVector(*atom["location"])
        if "charge" in atom:
            obatom.SetPartialCharge(atom["charge"])
    # If there is no bond data, try to infer them
    if "bonds" not in data or not data["bonds"]:
        obmol.ConnectTheDots()
        obmol.PerceiveBondOrders()
    # Otherwise, use the bonds in the data set
    else:
        for bond in data["bonds"]:
            if "atoms" not in bond:
                continue
            obmol.AddBond(bond["atoms"][0] + 1, bond["atoms"][1] + 1,
                          bond["order"])
    obmol.EndModify()
    if center:
        obmol.Center()
    return pybel.Molecule(obmol)


def pybel_to_json(molecule):
    """Converts a pybel molecule to json.

    Args:
        molecule: An instance of `pybel.Molecule`
    Returns:
        A Python dictionary containing atom and bond data
    """
    # Save atom element type and 3D location.
    atoms = [{"element": table.GetSymbol(atom.atomicnum),
              "location": atom.coords}
             for atom in molecule.atoms]
    # Saves partial charge data, if exists
    for atom, mol_atom in zip(atoms, molecule.atoms):
        if hasattr(mol_atom, "partialcharge") and mol_atom.partialcharge:
            atom["charge"] = mol_atom.partialcharge
    # Save number of bonds and indices of endpoint atoms
    bonds = [{"atoms": [b.GetBeginAtom().GetIndex(),
                        b.GetEndAtom().GetIndex()],
              "order": b.GetBondOrder()}
             for b in ob.OBMolBondIter(molecule.OBMol)]
    output = {"atoms": atoms, "bonds": bonds}

    # If there's unit cell data, save it to the json output
    if hasattr(molecule, "unitcell"):
        uc = molecule.unitcell
        output["periodic_connections"] = [[v.GetX(), v.GetY(), v.GetZ()]
                                          for v in uc.GetCellVectors()]
    return output


if __name__ == "__main__":
    # Lazy converter to test this out
    import sys
    in_data, in_format, out_format = sys.argv[1:]
    try:
        with open(in_data) as in_file:
            data = in_file.read()
    except IOError:
        data = in_data
    print convert(data, in_format, out_format, pretty=True)
