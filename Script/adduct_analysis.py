"""Metabolomics dataset analysis script.

Pipeline to find all instances of a metabolite according to its molecular mass and a 
given experimental M/Z ratio tolerance.

The script will use the cheatsheet.csv file to extract the rules used to compute the 
M/Z ratios of the potential adducts. If you want to use a different one, please keep
the same format as the one in this folder.

The datasets must be saved in the ../Datasets folder as .csv files. The names must 
include "neg" or "pos" somewhere, to let the script know whether the ions were
generated in positive or negative mode.

The output (if there are any matches) will be saved in the ../Output folder,
using as name <mass>_<tolerance>.csv.
The output indicates for each match the original file where it was found, the ionization
mode, the feature identifier, the retention time, the M/Z ratio, the rule used to 
compute the M/Z ratio of the adduct and whether the rule was used in direct or reverse
mode.

Steps:
    1. Load cheatsheet and extract rules
    2. For each dataset in ../Datasets:
        2.1. Compute M/Z of all potential adducts
        2.2. Compute upper and lower bounds (M/Z of adduct +/- tolerance)
        2.3. Check all M/Z in dataset and see if they fit in any ot the bounds
        2.4. Collect info on all features that fit in a bound
    3. Store info in a dataframe and save as .csv in ../Output
"""

import pandas as pd
import numpy as np
import os
from catcher import *
import argparse

############################################################################################

parser = argparse.ArgumentParser(description=__doc__,
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--mass',
                    default = 78.01393598,
                    type = float,
                    help="Metabolite mass")		

parser.add_argument('--tolerance',
                    default=0.002,
                    type = float,
                    help="M/Z tolerance to use when looking for adduct formation matches")

parser.add_argument("--cheatsheet",
                    default="cheatsheet.csv",
                    help="Name of the cheatsheet file to get the rules from")

args = parser.parse_args()

############################################################################################

def main(
        mass,
        tolerance,
        cheatsheet
        ):
    
    print("[analysis]: Run started")
    print(f"[analysis]: Using {cheatsheet} as adduct formation rules")
 
    cheatsheet = pd.read_csv(cheatsheet)
    ac = AdductCatcher()
    ac.get_rules(cheatsheet)
    
    filenames = os.listdir("../Datasets")
    paths = ["../Datasets/" + x for x in filenames]

    print(f"[analysis]: Found these files: {filenames}")

    match_file = []
    match_mode = []
    match_id = []
    match_rt = []
    match_mz = []
    match_rule = []
    match_dir = []
    
    for i in range(len(filenames)):
        print(f"[analysis]: Processing {filenames[i]}")
        db_target = pd.read_csv(paths[i])

        if "pos" in filenames[i]:
            mode = "Positive"
        else:
            mode = "Negative"

        temp_id, temp_rt, temp_mz, temp_rule, temp_dir = ac.process_db(mass,
                                                             tolerance,
                                                             db_target,
                                                             mode)

        if len(temp_id) > 0:            
            match_mode += [mode]*len(temp_id)
            match_file += [filenames[i]]*len(temp_id)
            match_id += temp_id
            match_rt += temp_rt
            match_mz += temp_mz
            match_rule += temp_rule
            match_dir += temp_dir
    
    print("[analysis]: Run finished")

    if len(match_file) > 0:
        output = pd.DataFrame({
            "File": match_file,
            "Mode": match_mode,
            "Feature ID": match_id,
            "Retention time": match_rt,
            "M/Z ratio": match_mz,
            "Adduct rule": match_rule,
            "Direction": match_dir
            })
        output_name = str(mass) + "_" + str(tolerance) + ".csv"
        output.to_csv("../Output/" + output_name)
        print(f"[analysis]: Run output saved as {output_name} in ../Output")
    else:
        print("[analysis]: No matches found, better luck next time :(")



if __name__ == "__main__":
    main(
        args.mass,
        args.tolerance,
        args.cheatsheet
        )







