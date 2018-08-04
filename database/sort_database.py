# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 23:26:20 2018

@author: Ahmad
"""
import os
import pandas as pd

os.chdir(r"C:\Users\Ahmad\Documents\GitHub\Chemical-Thermodynamics\database")
with open("inorganics.txt") as f:
    lines = f.readlines()

targetL = len(lines[0])
data = [lines[0].split(" ")]
data[0][-1] = data[0][-1][:-1] # remove \n from last element of the header
for line in lines[1:]:
    line = line.split(" ")
    L = len(line)
    
    for i,num in enumerate(line[-9:]):
        line[-(9 - i)] = float(num)
    
    while type(line[2]) == str:
        line = [line[0], line[1] + line[2]] + line[3:]
    data.append(line)

data = pd.DataFrame(data[1:],columns=data[0])
store = pd.HDFStore("database.h5")
store["inorganics"] = data
store.close()