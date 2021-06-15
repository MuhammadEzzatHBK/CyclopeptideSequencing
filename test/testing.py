# -*- coding: utf-8 -*-
"""
Created on Sun May 16 09:31:53 2021

@author: Muhammad Ayman Ezzat 
         Youmna Magdy Abdullah
"""
from algorithms import branch_and_bound
import timeit
import pandas as pd

''' Accuracy Testing '''
LabSpectrum = [97, 97, 99, 101, 103, 196, 198, 198, 200, 202, 295,
               297, 299, 299, 301, 394, 396, 398, 400, 400, 497]

LabResults = sorted(['PVCPT', 'PTPVC', 'PTPCV', 'PCVPT', 'VPTPC',
              'VCPTP', 'TPVCP', 'TPCVP', 'CPTPV', 'CVPTP'])

AssignmentResults = branch_and_bound(LabSpectrum)

print('Input: ', LabSpectrum)
print('Provided Lab Results: ', *LabResults)
print('Our Assignment Results: ', *AssignmentResults)
print('Are they identical? ', LabResults == AssignmentResults)

''' Perforamnce Testing '''
time_taken = []

for i in range(500):
    start = timeit.timeit()
    branch_and_bound(LabSpectrum)
    end = timeit.timeit()
    time_taken.append(abs(end - start))

data = {'duration' : time_taken}
DataFrame = pd.DataFrame(data)
DataFrame.to_csv('test_data.csv')