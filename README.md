# Adduct_finder
![Python 3.7](https://img.shields.io/badge/python-3.7%20%7C%203.8-brightgreen)  
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)  
Small script to find possible adducts of metabolites from a metabolomics MS analysis.  

## Dependencies:
- Numpy
- Pandas


## How does it work:
The script will first read possible adduct formation rules from the [cheatsheet.csv](./Script/cheatsheet.csv) file, depending on the ionization mode used in the MS analysis. Then, it will process each file in the [Datasets](Datasets) folder. First, it computes all possible M/Z ratios of the adducts associated to the given metabolite mass and ionization mode. Then, it computes upper and lower bounds of the M/Z ratios for each scenario according to the experimental tolerance. Then, it checks every M/Z ratio in the .csv files in the [Datasets](Datasets) folder to see if any fit within a boundary. If anything is a hit, it will collect the feature ID, from which table it came from, the retention time, the M/Z ratio, the ionization mode, the adduct formation rule and whether the rules were applied in direct or reverse mode. Once all hits are collected, they are saved in the [Output](Output) folder, as csv files following the metabolite mass + tolerance convention. Currently, there is an example output for DMSO using 0.002 tolerance.


## HOW TO:

To use the script, follow these steps:  

1. open a terminal instance (Linux) / anaconda prompt (Windows)   

2. Move to the Script folder (i.e. "cd Script")  

3. Write
``` 
python3 adduct_analysis.py --mass <metabolite_mass> --tolerance <tolerance> --cheatsheet <cheatsheet.csv>
```
   The words inside <...> are supposed to be changed according to your needs.
   For example, if I want to run an analysis for DMSO, using 0.002 as experimental tolerance and "example.csv" as
   adduct formation rules, here is the text I would need to write:
```
   python3 adduct_analysis.py --mass 78.01393598 --tolerance 0.002 --cheatsheet example.csv  
```
4. Check your outputs in the "Output" folder. If there were no matches, then there will be no file. The file will
   be named as <mass>_<tolerance>.csv

## NOTES:

- You can modify the cheatsheet file however you wish, as long as you write the coefficients, rules etc using the
same format for any new rules you wanna add. You can also remove rules if you want.  

- You can have multiple cheatsheet files in the scripts folder, just tell the script which one you wanna use for the
analysis using the "--cheatsheet" flag. Remember that the cheatsheet file has to be in the Script folder  

- The script will process automatically all files in the "Datasets" folder, make sure they are all .csv and have written
"pos" or "neg" somewhere so that the script knows which adduct formation rules to check.  

- The default values for --mass, --tolerance and --cheatsheet are 78.01393598, 0.002 and cheatsheet.csv. For example,
if you want to check a different metabolite but with same tolerance and rules, you can just write:  
```
python3 adduct_analysis.py --mass <your_new_mass>
```
If you don't specify the values, the script will automatically use the defaults.  

- If you want, you can also get instructions and a description of the steps the analysis does by writing in the terminal /
anaconda prompt the following line:
```
python3 adduct_analysis.py --help
```







