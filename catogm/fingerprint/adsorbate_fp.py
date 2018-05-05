# -*-coding: utf-8 -*-
# Author: Osman Mamun

import numpy as np
import json
import catogm


class Adsorbate_fp_generator():
    """class to generate adsorbate fingerprint for monoatomic species.
       This class is easily extensible to more complex molecules.
    """

    def __init__(self):
        """Initialize the class."""

    def _get_mendeleev_data(self):
        """Returns the Mendeleev data dictionary"""
        path = catogm.__file__.rsplit('/', 1)[0] + '/fingerprint/data/'
        return json.load(open(path + 'proxy-mendeleev.json'))

    def return_fp_list(self, atoms, fp_params):
        """TODO: documentation"""
        if not isinstance(atoms, object):
            raise NotImplementedError('{} data type not implemented.'.format(
                type(atoms)))

        fp_names = self.return_fp_names(fp_params)

        adsorbate_fp = []
        m_data = self._get_mendeleev_data()

        for k in fp_names:

            try:
                adsorbate_fp += [m_data[str(atoms.get_atomic_numbers()[0])]
                                       [k]]
            except TypeError:
                adsorbate_fp += [np.nan]

        return [i if i is not None else np.nan for i in adsorbate_fp]

    def return_fp_names(self, fp_params):
        """TODO: documentation"""
        fp_names = [i for i in fp_params]

        return fp_names
