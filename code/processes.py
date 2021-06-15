# -*- coding: utf-8 -*-
"""
Created on Sun May 16 09:31:53 2021

@author: Muhammad Ayman Ezzat 
         Youmna Magdy Abdullah
"""
from utils import get_molecular_weight, linear_spectrum, subset_spectrum, reversed_map
                    
def theoritical_spectrum(peptide_sequence):
   """Returns the theoritical spectrum of a given amino acid sequence.
    INPUT :
        peptide_sequence: string. The peptide sequence to get its theoritical spectrum
    OUTPUT:
        .: List. The theoritical spectrum of the given peptide sequence.
   """    
   linear_kmers = []
   cyclic_kmers = []
   for i in range(len(peptide_sequence)):
       for j in range(i,len(peptide_sequence)):
           linear_kmers.append(peptide_sequence[i:j+1])
   for i in range(2,len(peptide_sequence)):
       for j in range(i-1):
           cyclic_kmers.append(peptide_sequence[i:len(peptide_sequence)]+peptide_sequence[0:j+1])
   kmers =  linear_kmers+cyclic_kmers    
   return sorted(list(map(get_molecular_weight,kmers)))

def find_oneMers(spectrum):
   """Returns the one-mers that consist the sequence of the given spectrum.
    INPUT :
        spectrum: array-like. The spectrum required to get its one-mers.
    OUTPUT:
        .: List. A list of one-mers the consist the given spectrum.
    """    
   candidates = list('_'*len(spectrum))
   for i in range(len(spectrum)):
       if spectrum[i] in reversed_map:
           candidates[i] = reversed_map[spectrum[i]]
   return [cantdidate for cantdidate in candidates if cantdidate != '_' ]

def isConsistant(spectrum,kmer):
   """Checks whether a given kmer is consistent with a given spectrum or not.
    INPUT :
        spectrum: array-like. The spectrum required to check the given kmer against.
        kmer: string. The given kmer required to check its consistency.
    OUTPUT:
        .: bool. The consistency of the kmer against the spectrum.
    """     
   return subset_spectrum(spectrum,linear_spectrum(kmer))
    
def extend(oneMers,kmer):
   """Extends a given kmer.
    INPUT :
        oneMers: array-like. A list of amino acids that will be used to extend the given kmer.
        kmer: string. The given kmer required to extend.
    OUTPUT:
        .: List. A list of extended kmers.
    """        
   extentions = []
   for oneMer in oneMers:
       if(oneMers.count(oneMer) > list(kmer).count(oneMer)):
           extentions.append(kmer+oneMer)
   return extentions

