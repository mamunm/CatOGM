#!/usr/bin/env python
# -*-coding: utf-8 -*-

#bulk_fp_test.py
#Osman Mamun
#LAST UPDATED: 05-03-2018

from catogm.fingerprint.slab_fp import Slab_fp_generator
from ase.build import fcc111

slab = fcc111('Al', size=(2,2,3), vacuum=10.0)
slab1 = fcc111('Pd', size=(2,2,3), vacuum=10.0)
slab2 = fcc111('Pt', size=(2,2,3), vacuum=10.0)

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


sfp_gen = Slab_fp_generator()

list_slabs = [slab, slab1, slab2]
print(sfp_gen.return_fp(list_slabs, convoluted_params))
print(sfp_gen.return_fp_names(convoluted_params))

