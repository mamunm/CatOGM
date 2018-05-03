#!/usr/bin/env python
# -*-coding: utf-8 -*-

#bulk_fp_test.py
#Osman Mamun
#LAST UPDATED: 05-03-2018

from catogm.fingerprint.bulk_fp import Bulk_fp_generator
from ase.build import bulk

atoms = bulk('Pd', cubic=True)
atoms.set_chemical_symbols('Pd3Pt')

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

print(bfp_gen.return_fp_vec(atoms, convoluted_params, nonc_params, node=2))

