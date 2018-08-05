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
        self.molfracs = molfracs
        self.T = T # K
        self.P = P # bar
        
        # get property data and 
        self.data = self.search_db(self.compounds)
        self.params = self.mix_rules()
        self.v = self.pengrob_v()
    
    def set_state(self,T,P):
        self.T = T
        self.P = P
        self.params = self.mix_rules()
        self.v = self.pengrob_v()
        
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
        comps = store["all"]
        os.chdir(exec_path)
        store.close()
        
        df = pd.DataFrame()
        for name in names:
            comp = comps[comps["Name"].str.match(name,case=False)]
            if comp.shape[0] == 0:
                warnings.warn("No match found for component: %s" % name)
            elif comp.shape[0] > 1:
                warnings.warn("More than 1 match found for component: %s" % name)
            df = df.append(comp)
        
        return df
    
    def mix_rules(self):
        # get critical properties
        crits = self.data[["Tc[K]","Pc[bar]","omega"]]
        
        # a_i, b_i, and omega_i vectors
        aalpha_i = []; b_i = [];
        for _,Tc,Pc,omega in crits.itertuples():
            a_i = 0.45724*(R**2)*(Tc**2)/Pc
            b_i.append(0.07780*R*Tc/Pc)
            
            kappa = 0.37464 + 1.54226*omega - 0.26992*(omega**2)
            Tr = self.T/Tc
            alpha_i = (1 + kappa*(1 - np.sqrt(Tr)))**2
            aalpha_i.append(a_i*alpha_i)
        
        # aalpha_ii array
        aalpha_ii = []
        for aa_i in aalpha_i:
            row = []
            for aa_j in aalpha_i:
                row.append(np.sqrt(aa_i * aa_j))
            aalpha_ii.append(row)
        
        # sum according to mixing rules
        aalpha_mix = 0
        for i,y_i in enumerate(self.molfracs):
            for j,y_j in enumerate(self.molfracs):
                aalpha_mix += y_i*y_j*aalpha_ii[i][j]
        b_mix = sum([yi*bi for yi,bi in zip(self.molfracs,b_i)])
        
        return aalpha_mix,b_mix
    
    def pengrob_v(self):
        aalpha,b = self.params
        
        a3 = self.P
        a2 = self.P*b - R*self.T
        a1 = aalpha - 3*self.P*(b**2) - 2*R*self.T*b
        a0 = self.P*(b**3) + R*self.T*(b**2) - aalpha*b
        
        sols = np.roots([a3,a2,a1,a0])
        sols = np.real(sols[np.isreal(sols)]) # take only the real roots
        if len(sols) == 1:
            return [sols[0]]
        elif len(sols) != 0:
            print("Liquid-Vapour Region Detected!")
            return [min(sols), max(sols)]
        else:
            warnings.warn("Warning! No solution found for specific volume.")
            return
        
    
    
compounds = ["ethane","methane"]
mixture = fluid(compounds, [0.5,0.5], 273, 10)