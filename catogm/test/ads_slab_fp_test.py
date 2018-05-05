#!/usr/bin/env python
# -*-coding: utf-8 -*-

#bulk_fp_test.py
#Osman Mamun
#LAST UPDATED: 05-03-2018

from catogm.fingerprint.adsorbate_slab_fp import Adsorbate_slab_fp_generator
from ase.build import fcc111
from ase.build import add_adsorbate

slab = fcc111('Pd', size=(2,2,3), vacuum=10.0)
add_adsorbate(slab, 'C', 1.5, 'ontop')

ads_metal_params = ['atomic_number',                                            
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
                     'ads_connectivity']                                               
                                
metal_params = ['d-band_center',                                            
                'd-band_width',                                             
                'd-band_skewness',                                          
                'd-band_kurtosis']     


asfp_gen = Adsorbate_slab_fp_generator()


print(asfp_gen.return_fp_list(slab, ads_metal_params, metal_params))
print(asfp_gen.return_fp_names(ads_metal_params, metal_params))

