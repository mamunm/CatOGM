#!/usr/bin/env python
# -*-coding: utf-8 -*-

#bulk_fp_test.py
#Osman Mamun
#LAST UPDATED: 05-03-2018

from catogm.fingerprint.adsorbate_fp import Adsorbate_fp_generator
from ase.atoms import Atoms

atoms1 = Atoms('C')

fp_params = ['atomic_number',                                            
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
             'vdw_radius']                                            


afp_gen = Adsorbate_fp_generator()

print(afp_gen.return_fp_list(atoms1, fp_params))
print(afp_gen.return_fp_names(fp_params))

