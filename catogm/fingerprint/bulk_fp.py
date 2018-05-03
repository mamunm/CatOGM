#!/usr/bin/env python
# -*-coding: utf-8 -*-

#bulk_fp.py
#Osman Mamun
#LAST UPDATED: 05-02-2018

import numpy as np
import json
from ase.data import atomic_numbers
import pandas as pd
from catgen.utils import get_voronoi_neighbors

class Bulk_fp_generator():
    """class to generate bulk fingerprint for A1, L10, and L12 structures"""
    
    def __init__(self, atoms, convoluted_params, nonc_params, node=2):
        """Initialize the class.

        Parameters
        ----------
        atoms : Atoms object
              The atoms object to find fingerprints for.
        """
        if not isinstance(atoms, object):
            raise NotImplementedError('{} data type not implemented.'.format(
                type(atoms)))

        self.node = node
        if self.node > 2:
            raise NotImplementedError('node greater 2 is not implemented.')

    def _get_connectivity(self, atoms):
        """Returns the connectivity matrix. catgen.utils.get_voronoi_neighbors
        is used to get the connctivity matrix
        """
        return get_voronoi_neighbors(atoms)

    def _get_bulk_d_data(self):
        """Returns the bulk d band data dictionary"""

        return np.load('bulk_d_band_data.npy', encoding='latin1')[()]

    def _get_mendeleev_data(self):
        """Returns the Mendeleev data dictionary"""

        return json.load(open('proxy-mendeleev.json'))

    def _get_node(self, atoms):
        """Returns the node dictionary"""

        con = self._getconnectivity(atoms)
        node_dict = {}
        node_dict[0] = [a.symbol for a in atoms]

        node_dict[1] = []
        
        for i, c in enumerate(con):
            for j, cc in enumerate(c):
                if cc != 0:
                    node_dict[1].extend(cc * [[atoms[i].symbol, 
                                               atoms[j].symbol]])

        return node_dict

    def return_df(self, atoms, convoluted_params, nonc_params):

        """
        features = ['stoichiometry', 
                    'atomic_number', 
                    'atomic_radius', 
                    'atomic_volume',
                    'atomic_weight',
                    'boiling_point',
                    'covalent_radius_cordero',
                    'dipole_polarizability',
                    'electron_affinity',
                    'en_pauling',
                    'en_allen',
                    'en_ghosh',
                    'evaporation_heat',
                    'fusion_heat',
                    'group_id',
                    'period',
                    'heat_of_formation',
                    'lattice_constant_a',
                    'lattice_constant_c',
                    'melting_point',
                    'metallic_radius',
                    'specific_heat',
                    'thermal_conductivity',
                    'vdw_radius',
                    'd-band_center',
                    'd-band_width',
                    'd-band_skewness',
                    'd-band_kurtosis']
        """
        fp_list = [i for i in nonc_params]
        fp_list += [i + '_' + j for j in range(2) for i in convoluted_params]


   

