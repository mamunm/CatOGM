#!/usr/bin/env python
# -*-coding: utf-8 -*-

#bulk_fp.py
#Osman Mamun
#LAST UPDATED: 05-02-2018

import os
import numpy as np
import json
from ase.data import atomic_numbers
import pandas as pd
from catgen.utils import get_voronoi_neighbors
from ase.data import atomic_numbers as an
import multiprocessing
from tqdm import tqdm

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
        
        import catogm
        path = catogm.__file__.rsplit('/', 1)[0] + '/fingerprint/data/'
        return np.load(path + 
                       'bulk_d_band_data.npy', encoding='latin1')[()]

    def _get_mendeleev_data(self):
        """Returns the Mendeleev data dictionary"""

        import catogm
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

    def return_fp_list(self, atoms, convoluted_params, nonc_params):
        
        if not isinstance(atoms, object):
            raise NotImplementedError('{} data type not implemented.'.format(
                type(atoms)))
        

        fp_names= self.return_fp_names(convoluted_params, 
                                       nonc_params, io_mode='dict')

        bulk_fp = []
        
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
                bulk_fp += [atoms.cell[0][0]]
            
            if k == 'lattice_constant_c':
                bulk_fp += [atoms.cell[2][2]]

        
        node_dict = self._get_node(atoms)
        m_data = self._get_mendeleev_data()
        d_data = self._get_bulk_d_data()

        for k in fp_names['convoluted_params']:
            
            if all(['0' in k, not 'd-band' in k]):
                try:
                    bulk_fp += [sum([m_data[str(an[i])][k.rsplit('_', 1)[0]]**2 
                                              for i in node_dict[0]])]
                except TypeError:
                    bulk_fp += [np.nan]
            
            if all(['0' in k, 'd-band' in k]):
                try:
                    bulk_fp += [sum([d_data[i][0][k.split('_')[1]]**2 
                                              for i in node_dict[0]])]
                except (TypeError, KeyError) as e:
                    bulk_fp += [np.nan]

            if all(['1' in k, not 'd-band' in k]):
                try:
                    bulk_fp += [sum([m_data[str(an[i[0]])][k.rsplit('_', 1)[0]] 
                            * m_data[str(an[i[1]])][k.rsplit('_', 1)[0]] 
                            for i in node_dict[1]])]
                except TypeError:
                    bulk_fp += [np.nan]

            if all(['1' in k, 'd-band' in k]):
                try:
                    bulk_fp += [sum([d_data[i[0]][0][k.split('_')[1]] *
                                     d_data[i[1]][0][k.split('_')[1]]
                                              for i in node_dict[1]])]
                except (TypeError, KeyError) as e:
                    bulk_fp += [np.nan]
        
        return bulk_fp

    def return_fp_names(self, convoluted_params, 
                              nonc_params, io_mode='dict'):

        io_mo = ['dict', 'list']

        if io_mode not in io_mo:
            msg = 'Only dict and list type output mode is allowed.'
            raise NotImplementedError(msg)
        
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
   
    def return_fp(self, list_atoms, convoluted_params, nonc_params):
        
        if not isinstance(list_atoms, (list, object)):
            raise NotImplementedError('{} data type not implemented.'.format(
                type(list_atoms)))
        
        if not isinstance(list_atoms, list):
            return np.asarray(self.return_fp_list(list_atoms, 
                                                   convoluted_params, 
                                                   nonc_params))
        
        # Check for parallelized feature generation.
        fp_vector = []

        args = tuple((atoms, convoluted_params, nonc_params)
                     for atoms in list_atoms)

        if self.nprocs != 1:

            pool = multiprocessing.Pool(self.nprocs)
            parallel_iterate = pool.map_async(
                self.return_fp_list, args, callback=fp_vector.append)
            parallel_iterate.wait()
        
        else:
            
            for atoms in tqdm(list_atoms):
                fp_vector.append(self.return_fp_list(atoms,                  
                                                      convoluted_params,           
                                                      nonc_params))

        return np.asarray(fp_vector)
