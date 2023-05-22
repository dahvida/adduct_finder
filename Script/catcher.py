import numpy as np
import pandas as pd
from typing import *

class AdductCatcher:
    
    """
    Implementation of a class to catch all relevant features
    in a metabolomics dataset, given a mass and experimental tolerance.

    Contains all methods necessary to:
    1. Read all adduct formation rules from a cheatsheet file
    2. Compute all possible adducts given a molecular mass
    3. Check if any of the M/Z ratios measured in a run fit any
       of the possible adducts, according to a given experimental
       tolerance

    Attributes:
        coeffs:         Coefficients to use to convert metabolite mass
                        to M/Z of adduct
        shifts:         Mass shifts to use to convert metabolite mass
                        to M/Z of adduct
        names:          Rule for adduct formation
        modes:          Whether the rule applies for positive or negative
                        ion mode
        direction:      Whether the rule is used in direct or reverse mode
        n_rules:        Number of rules from the cheatsheet
        idx_pos:        Position of valid rules for positive mode
        idx_neg:        Position of valid rules for negative mode
    """

    def __init__(self):
        print("[ac]: Adduct Catcher created")
    
    def get_rules(self,
                  db: pd.DataFrame):

        #extract relevant columns from dataframe
        db.dropna(inplace=True)
        n_rules = len(db)
        names_1 = list(db["Ion name"])
        direction_1 = ["Direct"]*len(names_1)
        direction_2 = ["Reverse"]*len(names_1)
        modes = list(db["Mode"])
        charges = list(db["Mult"])
        shifts = list(db["Mass"]) 

        #prealloc coeffs and shifts
        self.coeffs = np.zeros((n_rules))
        self.shifts = np.zeros((n_rules))
        
        #get numerically accurate coeffs and shifts
        for i in range(n_rules):
            self.coeffs[i] = np.abs(np.float(charges[i]))
            if self.coeffs[i] == 0.33:
                self.coeffs[i] = 0.333333333333333333333333333333333333
            self.shifts[i] = np.float(shifts[i])  
        
        #adjust for dealing with direct and reverse mode
        self.coeffs = np.concatenate((self.coeffs, self.coeffs), axis=0)
        self.shifts = np.concatenate((self.shifts, -self.shifts), axis=0)
        self.names = names_1 + names_1
        self.modes = modes + modes
        self.direction = direction_1 + direction_2
        self.n_rules = len(self.names)
        
        #find indexes for positive/negative rules
        self.idx_pos = []
        self.idx_neg = []
        for i in range(len(self.modes)):
            if self.modes[i] == "Positive":
                self.idx_pos.append(i)
            else:
                self.idx_neg.append(i)
        
        print("[ac]: Rules set for the Adduct Catcher")
    
    
    def get_value(self,
                  mass: np.float,
                  rule_n: int):

        #get adduct M/Z
        return mass * self.coeffs[rule_n] + self.shifts[rule_n]
    
    
    def get_possibilities(self, mass, mode):
        #select right rules
        if mode == "Positive":
            idx = self.idx_pos
        else:
            idx = self.idx_neg
        
        #compute all M/Z according to ruleset
        results = np.zeros((len(idx)))
        for i in range(len(results)):
            results[i] = self.get_value(mass, idx[i])
        
        return results
    
    
    def get_match(self, mz, upper, lower):
        
        #check if M/Z is within any upper/lower range
        verdict = [False]*len(upper)
        for i in range(len(verdict)):
            if mz > lower[i]:
                if mz < upper[i]:
                    verdict[i] = True
        
        return verdict

        
    def process_db(self, mass, tolerance, db, mode):

        #get relevant info from target db
        rt = np.array(db)[:,2]
        mz = np.array(db)[:,3]
        names = list(db['General.All.ID'])
        
        #get relevant rules
        if mode == "Positive":
            mode_index = self.idx_pos
        else:
            mode_index = self.idx_neg
        
        #get M/Z boundaries for all rules
        upper = self.get_possibilities(mass, mode) + tolerance
        lower = self.get_possibilities(mass, mode) - tolerance
        
        #create output containers
        match_id = []
        match_rt = []
        match_mz = []
        match_rule = []
        match_direction = []

        #loop over all m/z
        for i in range(len(mz)):
            #get matches and their position
            verdicts = self.get_match(mz[i], upper, lower)
            match_index = np.where(verdicts)[0]

            #if there are matches, store info
            if len(match_index) > 0:
                for j in range(len(match_index)):
                    match_id.append(names[i])
                    match_rt.append(rt[i])
                    match_mz.append(mz[i])
                    match_rule.append(self.names[mode_index[match_index[j]]])
                    match_direction.append(self.direction[mode_index[match_index[j]]])

        return match_id, match_rt, match_mz, match_rule, match_direction
