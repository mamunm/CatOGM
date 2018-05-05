# -*-coding: utf-8 -*-
# Author: Osman Mamun

import numpy as np
import json
from catgen.utils import get_voronoi_neighbors
from ase.data import atomic_numbers as an
from tqdm import tqdm
import catogm


class Slab_fp_generator():
    """class to generate slab fingerprint for A1, L10, and L12 structures"""

    def __init__(self):
        """Initialize the class."""

    def _get_connectivity(self, atoms):
        """Returns the connectivity matrix."""
        return get_voronoi_neighbors(atoms)

    def _get_slab_d_data(self):
        """Returns the slab d band data dictionary"""
        path = catogm.__file__.rsplit('/', 1)[0] + '/fingerprint/data/'
        return np.load(path +
                       'slab_d_band_data.npy', encoding='latin1')[()]

    def _get_mendeleev_data(self):
        """Returns the Mendeleev data dictionary"""
        path = catogm.__file__.rsplit('/', 1)[0] + '/fingerprint/data/'
        return json.load(open(path + 'proxy-mendeleev.json'))

    def _get_node(self, atoms):
        """Returns the node dictionary"""
        con = self._get_connectivity(atoms)
        node_dict = {}
        node_dict[0] = [a.symbol for a in atoms]

        node_dict[1] = []

        for i, c in enumerate(con):
            for j, cc in enumerate(c):
                if cc != 0:
                    node_dict[1].extend(cc * [[atoms[i].symbol,
                                               atoms[j].symbol]])
        return node_dict

    def return_fp_list(self, atoms, convoluted_params):
        """TODO: documentation"""
        if not isinstance(atoms, object):
            raise NotImplementedError('{} data type not implemented.'.format(
                type(atoms)))

        fp_names = self.return_fp_names(convoluted_params)

        slab_fp = []

        node_dict = self._get_node(atoms)
        m_data = self._get_mendeleev_data()
        d_data = self._get_slab_d_data()

        for k in fp_names:

            if all(['0' in k, 'd-band' not in k]):
                try:
                    slab_fp += [sum([m_data[str(an[i])][k.rsplit('_', 1)[0]]**2
                                     for i in node_dict[0]])]
                except TypeError:
                    slab_fp += [np.nan]

            if all(['0' in k, 'd-band' in k]):
                try:
                    slab_fp += [sum([d_data[i][0][k.split('_')[1]]**2
                                     for i in node_dict[0]])]
                except (TypeError, KeyError) as e:
                    slab_fp += [np.nan]

            if all(['1' in k, 'd-band' not in k]):
                try:
                    slab_fp += [sum([m_data[str(an[i[0]])][k.rsplit('_', 1)[0]] *
                                     m_data[str(an[i[1]])][k.rsplit('_', 1)[0]]
                                     for i in node_dict[1]])]
                except TypeError:
                    slab_fp += [np.nan]

            if all(['1' in k, 'd-band' in k]):
                try:
                    slab_fp += [sum([d_data[i[0]][0][k.split('_')[1]] *
                                     d_data[i[1]][0][k.split('_')[1]]
                                     for i in node_dict[1]])]
                except (TypeError, KeyError) as e:
                    slab_fp += [np.nan]

        return slab_fp

    def return_fp_names(self, convoluted_params):
        """TODO: documentation"""
        fp_names = [i + '_' + str(j)
                    for i in convoluted_params
                    for j in range(2)]

        return fp_names

    def return_fp(self, list_atoms, convoluted_params):
        """TODO: documentation"""
        if not isinstance(list_atoms, (list, object)):
            raise NotImplementedError('{} data type not implemented.'.format(
                type(list_atoms)))

        if not isinstance(list_atoms, list):
            return np.asarray(self.return_fp_list(list_atoms,
                                                  convoluted_params))

        # Check for parallelized feature generation.
        fp_vector = []

        for atoms in tqdm(list_atoms):
            fp_vector.append(self.return_fp_list(atoms,
                                                 convoluted_params))

        return np.asarray(fp_vector)
