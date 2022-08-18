"""
pedigree_manipulator.py
> takes in the original all_calls.ped file
> pulls out each PARTID and finds the corresponding LI, RI, and SRS data from sample_summary.tsv
> outputs a pedigree file with ASD, LI, RI, and SRS genotype information
"""

# METHOD TO ADD LI, RI, AND SRS DATA TO PEDIGREE CHART
def add_disorders(ped_input, summary_input, ped_output):
    print("add_disorders method running...")
    
    # import modules
    import pandas as pd
    
    # import tab delimited files as pandas dataframes
    df_oldped = pd.read_csv(ped_input, sep = "\t")
    df_summary = pd.read_csv(summary_input, sep = "\t")
    
    # use PARTIDs from df_oldped to find corresponding data in df_summary
    df_oldped.rename(columns={"Phenotype": "ASD", "Individual ID": "PARTID"}, inplace = True)
    
    df_newped = df_oldped.merge(df_summary[["PARTID","LI", "RI", "SRS Cat"]], on="PARTID", how="left")
    df_newped.rename(columns={"SRS Cat":"SRS"}, inplace = True)
    
    # output merged dataframe as csv
    df_newped.to_csv(ped_output, sep = "\t")
    
    print("add_disorders method has completed running.")
    
# METHOD TO GENERATE SPECIFIC ASD, LI*, AND RI* PEDIGREES
def gen_spec_ped(ped_input, ASD_ped, ASD_LI_ped, ASD_RI_ped):
    print("gen_spec_ped() method is working...")
    
    # import modules
    import pandas as pd
    import numpy as np
    
    # import tab delimited files as pandas dataframes
    df_oldped = pd.read_csv(ped_input, sep = "\t", dtype = object, index_col = False)
    
    # generate ASD only pedigree
    df_ASD_only = df_oldped.iloc[:, 0:6]
    df_ASD_only = df_ASD_only.rename(columns = {"ASD":"Phenotype"})
    df_ASD_only.to_csv(ASD_ped, sep = "\t", index = False)
    print("ASD_only.ped file generated.")
    
    # generate LI only pedigree
    ped_file = open(ped_input, "r")
    ASD_LI_out = open(ASD_LI_ped, "w")
    ped_list = ped_file.read().split("\n")
    header_list = ped_list[0].split("\t")
    header_list[5] = "Phenotype"
    ASD_LI_out.write("\t".join(header_list[0:6]) + "\n")
    LI_star_counter = 0
    for a in range(1, len(ped_list)-2):
        line_list = ped_list[a].split("\t")
        if(line_list[6] == "2"):
            LI_star_counter += 1
            line_list[5] = "2"
        else:
            pass
        ASD_LI_out_string = "\t".join(line_list[0:6])
        ASD_LI_out.write(ASD_LI_out_string + "\n")
    ped_file.close()
    ASD_LI_out.close()
    print("ASD_LI.ped file generated. LI* individuals found: " + str(LI_star_counter))
    
    # generate RI only pedigree
    ped_file = open(ped_input, "r")
    ASD_RI_out = open(ASD_RI_ped, "w")
    ped_list = ped_file.read().split("\n")
    header_list = ped_list[0].split("\t")
    header_list[5] = "Phenotype"
    ASD_RI_out.write("\t".join(header_list[0:6]) + "\n")
    RI_star_counter = 0
    for a in range(1, len(ped_list)-2):
        line_list = ped_list[a].split("\t")
        if(line_list[7] == "2"):
            RI_star_counter += 1
            line_list[5] = "2"
        else:
            pass
        ASD_RI_out_string = "\t".join(line_list[0:6])
        ASD_RI_out.write(ASD_RI_out_string + "\n")
    ped_file.close()
    ASD_RI_out.close()
    print("ASD_RI.ped file generated. RI* individuals found: " + str(RI_star_counter))

    
    print("gen_spec_ped() method has finished working.")


# COMMANDS
# add_disorders("../Data/all_calls.ped", "../Data/sample_summary.tsv", "../Data/NJLAGS_CNV.ped")
gen_spec_ped("../Data/NJLAGS_CNV.ped", "../Data/ASD_only.ped", "../Data/ASD_LI.ped", "../Data/ASD_RI.ped")






