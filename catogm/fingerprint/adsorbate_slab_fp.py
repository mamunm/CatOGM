# -*-coding: utf-8 -*-
# Author: Osman Mamun

import numpy as np
import json
from catgen.utils import get_cutoff_neighbors
import pkg_resources


class Adsorbate_slab_fp_generator():
    """Class to generate adsorbate-slab fingerprint for monoatomic species.
    This class is easily extensible to more complex molecules.
    """

    def __init__(self):
        """Initialize the class."""

    def _get_connectivity(self, atoms):
        """Returns the connectivity matrix."""
        return get_cutoff_neighbors(atoms)

    def _get_mendeleev_data(self):
        """Returns the Mendeleev data dictionary."""
        data_path = pkg_resources.resource_filename(
            'catogm', 'data/proxy-mendeleev.json')
        print(data_path)

        with open(data_path) as f:
            data = json.load(f)

        return data

    def _get_slab_d_data(self):
        """Returns the slab d band data dictionary"""
        path = catogm.__file__.rsplit('/', 1)[0] + '/fingerprint/data/'
        return np.load(path +
                       'slab_d_band_data.npy', encoding='latin1')[()]

    def return_fp_list(self, atoms, ads_metal_params, metal_params):
        """TODO: documentation"""
        if not isinstance(atoms, object):
            raise NotImplementedError('{} data type not implemented.'.format(
                type(atoms)))

        fp_names = self.return_fp_names(ads_metal_params, metal_params)
        con = self._get_connectivity(atoms)
        con = con[-1][:-1]
        ads_an = atoms.get_atomic_numbers()[-1]
        slab_an = atoms.get_atomic_numbers()[:-1]
        slab_cs = atoms.get_chemical_symbols()[:-1]

        adsorbate_slab_fp = []

        d_data = self._get_slab_d_data()
        m_data = self._get_mendeleev_data()
        ads_con = np.count_nonzero(con)

        for k in fp_names:

            if all(['d-band' not in k, 'connectivity' not in k]):
                try:
                    adsorbate_slab_fp += [sum([c * m_data[str(slab_an[i])][k] *
                                               m_data[str(ads_an)][k]
                                               for i, c in enumerate(con)]) / ads_con]
                except TypeError:
                    adsorbate_slab_fp += [np.nan]

            if 'connectivity' in k:
                adsorbate_slab_fp += [ads_con]

            if 'd-band' in k:
                try:
                    adsorbate_slab_fp += [sum([c * d_data[slab_cs[i]][0]
                                               [k.split('_')[1]]
                                               for i, c in enumerate(con)]) / ads_con]
                except TypeError:
                    adsorbate_slab_fp += [np.nan]

        return adsorbate_slab_fp

    def return_fp_names(self, ads_metal_params, metal_params):
        """TODO: documentation"""
        fp_names = [i for i in ads_metal_params]
        fp_names += [i for i in metal_params]

        return fp_names
