#!/usr/bin/env python
# -*-coding: utf-8 -*-

#bulk_fp_test.py
#Osman Mamun
#LAST UPDATED: 05-03-2018

from catogm.fingerprint.bulk_fp import Bulk_fp_generator
from ase.build import bulk

atoms1 = bulk('Pd', cubic=True)
atoms1.set_chemical_symbols('Al3Pt')
atoms2 = bulk('Pd', cubic=True)
atoms2.set_chemical_symbols('Pd3Pt')
atoms3 = bulk('Pd', cubic=True)

convoluted_params = features = ['atomic_number',                                            
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
                                'melting_point',                                            
                                'metallic_radius',                                          
                                'specific_heat',                                            
                                'thermal_conductivity',                                     
                                'vdw_radius',                                               
                                'd-band_center',                                            
                                'd-band_width',                                             
                                'd-band_skewness',                                          
                                'd-band_kurtosis']     

nonc_params = ['stoichiometry', 
               'lattice_constant_a',                                       
               'lattice_constant_c']                                       

bfp_gen = Bulk_fp_generator()

list_atoms = [atoms1, atoms2, atoms3]
print(bfp_gen.return_fp(list_atoms, convoluted_params, nonc_params))
print(bfp_gen.return_fp_names(convoluted_params, nonc_params, io_mode='list'))

