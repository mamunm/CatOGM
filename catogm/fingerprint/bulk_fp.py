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
    
    def __init__(self, nprocs=1):
        """Initialize the class.

        Parameters
        ----------
        nprocs : int
            Number of cores available for parallelization. Default is 1, e.g.
            serial. Set None to use all available cores.
        """
        self.nprocs = nprocs

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

    def return_fp_vec(self, atoms, convoluted_params, nonc_params, node):
        
        if not isinstance(atoms, object):
            raise NotImplementedError('{} data type not implemented.'.format(
                type(atoms)))

        fp_names= self.return_fp_names(convoluted_params, 
                                       nonc_params, node=2, io_mode='dict')

        fp_list = []
        
        for k in fp_names['nonc_params']:
            if k == 'stoichiometry':
                try: 
                    bulk_fp +=  [0.71 if atoms.info['SBS'] == 'L10'
                            else 0.79 if atoms.info['SBS'] == 'L12'
                            else 1.00]

                except KeyError:
                    chem_form = atoms.get_chemical_formula()
                    if '3' in chem_form:
                        bulk_fp += [0.79]
                    elif '2' in chem_form:
                        bulk_fp += [0.71]
                    else:
                        bulk_fp += [1.00]
            
            if k == 'lattice_constant_a':
                bulk_fp += atoms.cell[0][0]
            
            if k == 'lattice_constant_c':
                bulk_fp += atoms.cell[2][2]

        
        for k in fp_names['convoluted_params']:
            bulk_fp += _get_fp(k)
        



        return np.asarray(bulk_fp)

    def return_fp_names(self, convoluted_params, 
                              nonc_params, node=2, 
                              io_mode='dict'):

        io_mo = ['dict', 'list']

        if io.mode not in io_mo:
            raise NotImplementedError('Only dict and list type output mode
                    id allowed.')
        
        if node > 2:
            raise NotImplementedError('node greater 2 is not implemented.')

        if io_mode == 'dict':
            fp_names = {}
            fp_names['nonc_params'] = [i for i in nonc_params]
            fp_names['convoluted_params'] = [i + '_' + str(j) 
                                     for i in convoluted_params 
                                     for j in range(2)]

        if io_mode == 'list':
            fp_names = [i for i in nonc_params]
            fp_names += [i + '_' + str(j) 
                                     for i in convoluted_params 
                                     for j in range(2)]


        return fp_names

   

