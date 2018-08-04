# -*- coding: utf-8 -*-
"""
Created on Sat Aug  4 00:47:46 2018

@author: Ahmad
"""
import pandas as pd
import numpy as np
import os
import warnings

'''
Universal constants
'''
R = 8.314/100 # m^3*bar/(kmol*K)

class fluid(object):
    def __init__(self,compounds,molfracs,T,P):
        self.compounds = compounds
        self.T = T # K
        self.P = P # bar
        
        self.data = self.search_db(compounds)
    
    @staticmethod
    def search_db(names):
        '''
        This function searches the current database of compounds for:
            Molecular weights, criticial properties, acentricity, and antoine's
            equation parameters
        '''
        # retrieve compounds and parameters from database
        exec_path = os.getcwd()
        os.chdir("../database")
        store = pd.HDFStore("database.h5")
        comps = [store["organics"],store["inorganics"]]
        os.chdir(exec_path)
        store.close()
        comps = pd.concat(comps)
        comps.index = [i for i in range(0,len(comps))]
        
        df = pd.DataFrame()
        for name in names:
            comp = comps[comps["Name"].str.match(name,case=False)]
            if comp.shape[0] == 0:
                warnings.warn("No match found for component: %s" % name)
            elif comp.shape[0] > 1:
                warnings.warn("More than 1 match found for component: %s" % name)
            df = df.append(comp)
        
        return df
    
    
compounds = ["ethane","methane"]
mixture = fluid(compounds, [0.5,0.5], 273, 1)