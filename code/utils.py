# -*- coding: utf-8 -*-
"""
Created on Sun May 16 09:31:53 2021

@author: Muhammad Ayman Ezzat 
         Youmna Magdy Abdullah
"""
from itertools import chain

molecular_weights = {'G':57,'A':71,'S':87,
                    'P':97,'V':99,'T':101,
                    'C':103,'I':113,'L':113,
                    'N':114,'D':115,'K':128,
                    'Q':128,'E':129,'M':131,
                    'H':137,'F':147,'R':156,
                    'Y':163,'W':186}

reversed_map =  {v:k for k, v in molecular_weights.items()}

def get_molecular_weight(kmer):
   """Calculates molecular weight of a given amino acid sequence.
    INPUT :
        kmer: string. The amino acid sequence to get its weight
    OUTPUT:
        weight: integer. The molecular weight of the given amino acid sequence.
    """
   weight = 0
   for AA in kmer:
       weight += molecular_weights[AA]
   return weight

def linear_spectrum(peptide_sequence):
    """Returns the linear spectrum of a given amino acid sequence.
    INPUT :
        peptide_sequence: string. The peptide sequence to get its linear spectrum
    OUTPUT:
        .: List. The linear spectrum of the given peptide sequence.
    """    
    linear_kmers = []
    for i in range(len(peptide_sequence)):
        for j in range(i,len(peptide_sequence)):
            linear_kmers.append(peptide_sequence[i:j+1])
    return sorted(list(map(get_molecular_weight,linear_kmers)))

def subset_spectrum(spectrum,sub_spectrum):
   """Checks whether a sub spectrum is subset from a spectrum or not.
    INPUT :
        spectrum: array-like. The spectrum of the RHS.
        sub_spectrum: array-like. The spectrum of the LHS.
    OUTPUT:
        .: bool. Is sub_spectrum subset from spectrum.
   """    
   for weight in sub_spectrum:
       if weight not in spectrum:
           return False
   return True

def flatten(mapper_object):
   """Transforms a mapper object data layout.
    INPUT :
        mapper_object: mapper object. The mapper object required to transform.
    OUTPUT:
        .: List. The list transformation of the mapper object.
    """        
   return list(chain.from_iterable(list(mapper_object)))