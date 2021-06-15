# -*- coding: utf-8 -*-
"""
Created on Sun May 16 09:31:53 2021

@author: Muhammad Ayman Ezzat 
         Youmna Magdy Abdullah
"""
from functools import partial
from itertools import compress
from utils import flatten
from processes import theoritical_spectrum, find_oneMers, isConsistant, extend

def branch_and_bound(spectrum):
   """Our implementation of the Branch & Bound algorithm for the cyclopeptide sequencing problem.
    INPUT :
        spectrum: array-like. The input spectrum.
    OUTPUT:
        .: List. All the linear representations of the cyclic sequence of the protein (Peptide Sequence).
    """        
   oneMers = find_oneMers(spectrum)
   extend_kmers = partial(extend,oneMers)
   temp_list = flatten(map(extend_kmers,oneMers))
   isConsistant_prefilled = partial(isConsistant,spectrum)
   keepers = list(map(isConsistant_prefilled,temp_list))
   solutions = list(compress(temp_list,keepers))
   lengthes = [-1,-2]
   while(lengthes[-1] != lengthes[-2]):
       temp_list = flatten(map(extend_kmers,solutions))
       keepers = list(map(isConsistant_prefilled,temp_list))
       solutions = list(set(list(compress(temp_list,keepers))))
       lengthes.append(len(solutions))
   theoritical_spectra = list(map(theoritical_spectrum,solutions))
   validation = [theoritical_spectrum_ == spectrum for theoritical_spectrum_ in theoritical_spectra]
   return sorted(list(compress(solutions,validation)))


