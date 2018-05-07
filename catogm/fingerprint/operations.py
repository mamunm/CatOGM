import numpy as np


def raw_properties(
        atoms=None,
        atoms_properties=None,
        connectivity=None):
    """Return all atom properties without manipulation."""
    fingerprint = np.concatenate(atoms_properties)

    return fingerprint


def convolution_d0(
        atoms=None,
        atoms_properties=None,
        connectivity=None):
    """Return the square of each property with each atom."""
    convolution = np.dot(atoms_properties,
                         atoms_properties.T).diagonal()

    return convolution


def convolution_d1(
        atoms=None,
        atoms_properties=None,
        connectivity=None):
    """Return the multiple of each property for all neighboring atom
    pairs.
    """
    V = np.dot(atoms_properties, connectivity)
    convolution = np.dot(V, atoms_properties.T).diagonal()

    return convolution


def bonding_convolution(
        atoms=None,
        atoms_properties=None,
        connectivity=None):
    """Perform convolution of metal atoms with bonded adsorbates."""
    # This is a CatKit convention
    bond_index = np.where(atoms.get_tags() == -1)[0]

    V = np.dot(atoms_properties, connectivity)[:, bond_index]
    P = np.dot(V, atoms_properties[:, bond_index].T).diagonal()
    convolution = P / connectivity[bond_index].sum()

    return convolution
