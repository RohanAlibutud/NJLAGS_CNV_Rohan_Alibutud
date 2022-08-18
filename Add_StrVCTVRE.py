"""
Add_StrVCTVRE.py
> takes an annotated .tsv file and extracts positions
> composes positions into .bed file
> inputs .bed file into StrVCTVRE 
"""

import pandas as pd
import numpy as np

print("Running Add_StrVCTVRE.py script...")

# ADD StrVCTVRE ANNO METHOD
def add_STR(STR_input, var_input, var_output):
    
    # open files
    STR_in = open(STR_input, "r")
    var_in = open(var_input, "r")
    var_out = open(var_output, "w")
    
    # read STR_in and var_in as pandas dataframes
    df_str = pd.read_csv(STR_in, sep = "\t", header = None)
    df_str.columns = ["chrom", "start", "stop", "svtype", "StrVCTVRE"]
    df_var = pd.read_csv(var_in, sep = "\t")
    df_var["StrVCTVRE"] = df_str["StrVCTVRE"].copy()
    print(df_var)
    
    df_var.to_csv(var_output, sep = "\t", index = False)
    
    # close files
    STR_in.close()
    var_in.close()
    var_out.close()
    
print("Add_StrVCTVRE.py has completed running.")

add_STR("../Data/GMN_StrVCTVRE.tsv", "../Data/GMNP_intermediate5.tsv", "../Data/GMNP_out_2022_03_01.tsv")