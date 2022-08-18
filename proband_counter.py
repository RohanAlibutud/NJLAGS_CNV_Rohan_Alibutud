"""
proband_counter.py
> Takes in tab delimited text file
> Excludes lines without genotype information
> Counts probands for ASD, LI, RI, SRS
"""

def count_probands(input_file):
    import pandas as pd
    
    print("count_probands method running...")
    
    df_summary = pd.read_csv(input_file, delimiter = "\t")
    df_proband = df_summary[(df_summary["VCF"] == 1)]
    print(df_proband["ASD"].value_counts())
    print(df_proband["ADHD_evidence"].value_counts())
    print(df_proband["LI"].value_counts())
    print(df_proband["RI"].value_counts())
    print(df_proband["SRS Cat"].value_counts())
    
    print("count_probands method has completed running.")
    
count_probands("sample_summary.tsv")