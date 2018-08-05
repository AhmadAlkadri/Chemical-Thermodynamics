# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 23:26:20 2018

@author: Ahmad
"""
import os
import pandas as pd

os.chdir(r"C:\Users\Ahmad\Documents\GitHub\Chemical-Thermodynamics\database")

def Critpropstxt_to_HDF5(fname):
    '''
    NOTE: update this to reflect column update before calling...    
    '''
    with open(fname) as f:
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
            line = [line[0], line[1] + " " + line[2]] + line[3:]
        data.append(line)
    
    data = pd.DataFrame(data[1:],columns=data[0])
    store = pd.HDFStore("database.h5")
    store["inorganics"] = data
    store.close()

def join_HDF5():
    store = pd.HDFStore("database.h5")
    comps = [store["organics"],store["inorganics"]]
    comps = pd.concat(comps)
    comps.index = [i for i in range(0,len(comps))]
    
    keep_cols = ['Formula', 'Name', 'MW[g/mol]', 'Tc[K]', 'Pc[bar]', 'omega']
    ant_cols = ['A', 'B', 'C', 'Tmin', 'Tmax']
    ant_params = [[A,B,C,Tmin,Tmax] for _,A,B,C,Tmin,Tmax in comps[ant_cols].itertuples()]
    ant_params = pd.Series(ant_params,name="Antoine Constants")
    new = pd.concat([comps[keep_cols], ant_params], axis=1)
    
    store["all"] = new
        
    store.close()

def Cpgastxt_to_HDF5(fname):
    '''
    Functional form is:
        Cp/R = A + B*T + C*T**2 + D*T**(-2) + E*T**3 with T [=] K
    '''
    with open(fname) as f:
        lines = f.readlines()
    
    targetL = len(lines[0])
    data = [lines[0].split(" ")]
    data[0][-1] = data[0][-1][:-1] # remove \n from last element of the header
    for line in lines[1:]:
        line = line.split(" ")
        L = len(line)
        
        for i,num in enumerate(line[-8:]):
            line[-(8 - i)] = float(num)
        
        while type(line[2]) == str:
            line = [line[0], line[1] + " " + line[2]] + line[3:]
        data.append(line)
    
    data = pd.DataFrame(data[1:],columns=data[0])
    store = pd.HDFStore("database.h5")
    store["/Cp/gases"] = data
    store.close()
    
Cpgastxt_to_HDF5("Cp_gas.txt")

Critpropstxt_to_HDF5("organics.txt")
Critpropstxt_to_HDF5("inorganics.txt")
join_HDF5()
    
    
    
    
    
    
    
    
    
    
    
