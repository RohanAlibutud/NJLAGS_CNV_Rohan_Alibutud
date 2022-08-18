"""
CNV_merge.py
> reformats merged .tsv files as a true vcf-like
> reads in all three files and compiles a list of individuals
> uses list of individuals to build a table with a separate column for each
> read through each merged .tsv in order to populate the table
"""

import os
import numpy as np

# METHOD TO COMPILE LIST OF INDIVIDUALS
def indv_find(input_folder):
    print("indv_find script running...")
    
    # open input files
    axiom_in = open(input_folder + "axiom_merged.tsv", "r")
    jan_in = open(input_folder + "jan_merged.tsv", "r")
    may_in = open(input_folder + "may_merged.tsv", "r")
    
    # open output file
    indv_out = open(input_folder + "individuals.txt", "w")
    
    # declare list for holding all the individual IDs
    indv_list = []
    
    # iterate through all the axiom individuals
    for line in axiom_in:
        line_list = line.split("\t")
        
        # skip header line
        if(line_list[3] == "partid"):
            continue
        
        indv_ID = int(line_list[3])
        # append only unique individual IDs
        if(indv_ID in indv_list):
            continue
        else:
            indv_list.append(indv_ID)
            
    # iterate through all the jan individuals
    for line in jan_in:
        line_list = line.split("\t")
        
        # skip header line
        if(line_list[3] == "partid"):
            continue
        
        indv_ID = int(line_list[3])
        # append only unique individual IDs
        if(indv_ID in indv_list):
            continue
        else:
            indv_list.append(indv_ID)
    
    # iterate through all the may individuals
    for line in may_in:
        line_list = line.split("\t")
        
        # skip header line
        if(line_list[3] == "partid"):
            continue
        
        indv_ID = int(line_list[3])
        # append only unique individual IDs
        if(indv_ID in indv_list):
            continue
        else:
            indv_list.append(indv_ID)
        
    # sort indv_list
    indv_list.sort()
    
    # write indv_list out to file
    for i in range(0,len(indv_list)):
        indv_list[i] = str(indv_list[i])
    indv_out.write("\n".join(indv_list))

    # close input files
    axiom_in.close()
    jan_in.close()
    may_in.close()
    
    # close output file
    indv_out.close()
    print("indv_find script complete.")

# METHOD TO CHECK IF CALLERS DIFFER IN COPY NUMBER
def check_CNV(input_file):
    file_in = open(input_file, "r")
    
    if("axiom" in input_file):
        filename = "axiom_merged_2"
    if("jan" in input_file):
        filename = "jan_merged_3"
    if("may" in input_file):
        filename = "may_merged_2"
    if("all" in input_file):
        filename = "all_merged_2"
    
    file_out = open("../Data/merge_outputs/" + filename + ".mismatches.txt", "w")
    
    # iterate through file line-by-line
    for line in file_in:
        
        # declare variable outside inline loop
        CNV = "ERROR"
        
        # separate line into tab-delimited list
        line_list = line.split("\t")
        
        # iterate through line column-by-column
        for i in range(0, len(line_list)):
            if("/" in line_list[i]): # upon finding a genotype filed
                if(CNV == "ERROR"): # replace placeholder with first CNV
                    CNV = line_list[i].split("/")[0]
                elif(CNV == line_list[i].split("/")[0]): # if matching, proceed as normal
                    continue
                else:
                    file_out.write(line) # if mismatched, write out to mismatch file
                    break
    
    file_in.close()
    
# METHOD TO SPLIT
def CNV_split(input_list, CNV_list, caller_mismatch):
    # separate deletion and duplication columns
    del_list = input_list.copy()
    dup_list = input_list.copy()
    
    # in a deletion variant:
    # recalculate outermost positions
    del_start = "0"
    del_stop = "0"
    for j in range(0, len(del_list)):
        # identify genotype column and delete duplications
        if("/" in del_list[j]):
            if(del_list[j].split("/")[0] == "2" or del_list[j].split("/")[0] == "3" or del_list[j].split("/")[0] == "4"):
                del_list[j] = "."
            else:
                # determine new del_start
                if(int(del_start) == 0):
                    del_start = del_list[j].split(":")[4]
                elif(int(del_start) != 0 and int(del_start) >= int(del_list[j].split(":")[4])):
                    del_start = del_list[j].split(":")[4]

                # determine new del_stop
                if(int(del_stop) == 0):
                    del_stop = del_list[j].split(":")[5]
                elif(int(del_stop) != 0 and int(del_stop) <= int(del_list[j].split(":")[5])):
                    del_stop = del_list[j].split(":")[5]
                
                
        # check if callers differ in CNV
        if("|" in del_list[j]):
            del_indv_list = del_list[j].split("|")
            qsnp_CNV = del_indv_list[0][0]
            penncnv_CNV = del_indv_list[1][0]
            if(qsnp_CNV != penncnv_CNV):
                caller_mismatch = True
            
            # assign avg del
            if(caller_mismatch == True):
                for k in range(0, len(del_indv_list)):
                    del_col_list = del_indv_list[k].split("/")
                    avg_del = round(((qsnp_CNV + penncnv_CNV)/2),1)
                    del_col_list[0] = str(avg_del)
                    del_indv_list[k] = "/".join(del_col_list)
            del_list[j] = "|".join(del_indv_list)
            
        # apply new del_start and del_stop
        del_list[1] = del_start
        del_list[2] = del_stop
    
    del_string = "\t".join(del_list)
                                        
    # in a duplication variant:
    # recalculate outermost positions
    dup_start = "0"
    dup_stop = "0"
    
    # calculate average CNV for duplications
    dup_CNVs = CNV_list.copy()
    dup_CNVs[:] = [x for x in dup_CNVs if x != "0"]
    dup_CNVs[:] = [x for x in dup_CNVs if x != "1"]
    dup_CNVs = [ int(x) for x in dup_CNVs ]
    
    
    for j in range(0, len(dup_list)):
        # identify genotype column and delete deletions
        if("/" in dup_list[j]):
            if(dup_list[j].split("/")[0] == "0" or dup_list[j].split("/")[0] == "1"):
                dup_list[j] = "."
            else:
                # determine new dup_start
                if(int(dup_start) == 0):
                    dup_start = dup_list[j].split(":")[4]
                elif(int(dup_start) != 0 and int(dup_start) >= int(dup_list[j].split(":")[4])):
                    dup_start = dup_list[j].split(":")[4]

                # determine new dup_stop
                if(int(dup_stop) == 0):
                    dup_stop = dup_list[j].split(":")[5]
                elif(int(dup_stop) != 0 and int(dup_stop) <= int(dup_list[j].split(":")[5])):
                    dup_stop = dup_list[j].split(":")[5]
            
            # set average CNV for each dup
           
            if("|" in dup_list[j]):
                dup_indv_list = dup_list[j].split("|")
                # check if callers differ in CNV
                qsnp_CNV = int(dup_indv_list[0][0])
                penncnv_CNV = int(dup_indv_list[1][0])
                if(qsnp_CNV != penncnv_CNV):
                    caller_mismatch = True
                
                # assign avg dup
                if(caller_mismatch == True):
                    for k in range(0, len(dup_indv_list)):
                        dup_col_list = dup_indv_list[k].split("/")
                        avg_dup = round(((qsnp_CNV + penncnv_CNV)/2),1)
                        dup_col_list[0] = str(avg_dup)
                        dup_indv_list[k] = "/".join(dup_col_list)
                dup_list[j] = "|".join(dup_indv_list)
                                              
                
        # apply new dup_start and dup_stop
        dup_list[1] = dup_start
        dup_list[2] = dup_stop
    
    # write out strings
    del_list[3] = "SPLIT, DELETIONS"
    del_string = "\t".join(del_list)
    dup_list[3] = "SPLIT, DUPLICATIONS"
    dup_string = "\t".join(dup_list)
    
    output_string_list = [del_string, dup_string]
    
    return(output_string_list)

# METHOD TO MERGE
def CNV_merge(input_file, output_file):
    print("CNV_merge script running...")
    
    # open input files
    file_in = open(input_file, "r")
    indv_in = open("../Data/individuals.txt", "r")
    
    # open output files
    output_filename = input_file.replace("_sorted.tsv","") + output_file + ".tsv"
    file_out = open(output_filename, "w")
    caller_mismatch_out = open("../Data/merge_outputs/caller_mismatches.tsv","w")
    
    error_out = open("error.txt", "w")
    
    # PRINT OUT HEADER
    # set individual list
    indv_list = indv_in.read().split("\n")
    header_list = ["#Chrom", "Start", "End", "SV_type"]
    header_list.extend(indv_list)
    header_list.append("\n")
    header_string = "\t".join(header_list)
    
    # write out header
    file_out.write(header_string)
    
    # declare ongoing merge variables
    info_list = []
    qsnp_index = 0
    penncnv_index = 0
    final_outer_start = 0
    final_outer_stop = 0
    merge_count = 0 # keeps track of how many merges are undergone
    
    # POPULATE TABLE
    # read in file as newline-separated list
    file_list = file_in.read().split("\n")
    
    # out-of-loop variables
    merge_ongoing = False # boolean that is used to enable multiple-line merges
    
    # run through the entire sorted_list, reading two lines at a time
    for i in range(1, len(file_list) - 1):
        
        # declare current line and next line as lists
        curr_line = file_list[i].split("\t")
        next_line = file_list[i + 1].split("\t")
        
        # skip any empty lines at the end
        if(len(next_line) <= 1 or len(curr_line) <= 1):
            continue
        
        # assign positions and determine length for current line
        curr_chrom = curr_line[0]
        curr_start = int(curr_line[1])
        curr_stop = int(curr_line[2])
        curr_ID = curr_line[3]
        curr_info = curr_ID + ":" + curr_line[4]
        curr_info_list = curr_info.split(":")
        curr_CNV = curr_info_list[5]
        curr_geno = [curr_CNV + "/.", ":".join(curr_info_list[0:5]), ":".join(curr_info_list[6:8])]
        curr_geno = ":".join(curr_geno)
        curr_length = curr_stop - curr_start
        
        # assign positions and determine length for next line
        next_chrom = next_line[0]
        next_start = int(next_line[1])
        next_stop = int(next_line[2])
        next_ID = next_line[3]
        next_info = next_ID + next_line[4]
        next_info_list = next_info.split(":")
        next_CNV = next_info_list[5]
        next_geno = [next_CNV + "/.", ":".join(next_info_list[0:5]), ":".join(next_info_list[6:8])]
        next_geno = ":".join(next_geno)
        next_length = next_stop - next_start
        
        
        # skip sex chromosomes
        if curr_chrom == "X" or curr_chrom == "Y":
            continue
        
        # declare variable to hold information of whether it's a mismatch and what kind
        split_string = "NO MISMATCH"
        
        # CHECK IF THEY OVERLAP SUCH THAT THEY SHOULD BE MERGED
        overlap = False
        
        
        
        # must first be from the same chromosome
        if(curr_chrom == next_chrom):
            
            # determine overlap (inner) and outer positions
            if(curr_start >= next_start): # OVERLAP START
                overlap_start = curr_start # greater start is rightmost or inner start
                outer_start = next_start # lesser start is leftmost or outer start
            else:
                overlap_start = next_start # greater start is rightmost or inner start
                outer_start = curr_start # lesser start is leftmost or outer start
                
            if(curr_stop <= next_stop): # OVERLAP STOP
                overlap_stop = curr_stop # lesser stop is the leftmost or inner stop
                outer_stop = next_stop # greater stop is rightmost or outer stop
            else:
                overlap_stop = next_stop # lesser stop is the leftmost or inner stop
                outer_stop = curr_stop # greater stop is the rightmost or outer stop
            
            # adjust final_outer variables
            if(final_outer_start == 0 or final_outer_start >= outer_start):
                final_outer_start = outer_start
            else:
                final_outer_start = final_outer_start
            if(final_outer_stop == 0 or final_outer_stop <= outer_stop):
                final_outer_stop = outer_stop
            else:
                final_outer_stop = final_outer_stop
            
            # calculate overlap length
            overlap_length = overlap_stop - overlap_start
            
            # perform overlap 70% calculation
            curr_70 = curr_length * .7
            next_70 = next_length * .7
            if(overlap_length >=  curr_70 and overlap_length >= next_70):
                overlap = True
            else:
                overlap = False
        
        # if chromosomes do not match, no overlap is possible
        else:
            overlap = False
            
        # if overlap is true, perform merge
        if(overlap == True):
            
            # PERFORM MERGE
            merge_ongoing = True # initiate merge process
            merge_count += 1
            
            # declare merge list for output
            merge_list = [curr_chrom, str(final_outer_start), str(final_outer_stop), split_string]
            
            # add individual columns
            for k in range(0, len(indv_list)):
                merge_list.append(".")
            
            # append caller information
            curr_geno_list = curr_geno.split(":")
            if("qsnp" in curr_geno_list[2]):
                qsnp_index += 1
                curr_geno_list[2] = curr_geno_list[2] + "_" + str(qsnp_index)
            elif("penncnv" in curr_geno_list[2]):
                penncnv_index += 1
                curr_geno_list[2] = curr_geno_list[2] + "_" + str(penncnv_index)
            curr_geno = ":".join(curr_geno_list)
            info_list.append(curr_geno)            
            
        else:
            # if there is an ongoing merge, output it first
            if(merge_ongoing == True):   
                
                # append caller information for last line in merge
                curr_geno_list = curr_geno.split(":")
                if("qsnp" in curr_geno_list[2]):
                    qsnp_index += 1
                    curr_geno_list[2] = curr_geno_list[2] + "_" + str(qsnp_index)
                elif("penncnv" in curr_geno_list[2]):
                    penncnv_index += 1
                    curr_geno_list[2] = curr_geno_list[2] + "_" + str(penncnv_index)
                curr_geno = ":".join(curr_geno_list)
                info_list.append(curr_geno)
                
                # add info_list from current individual to the corresponding column
                # iterate through each call at this position
                for a in range(0, len(info_list)):
                    # isolate ID
                    info_ID = info_list[a].split(":")[1]
                    # add caller information to the corresponding individual column
                    for b in range(0, len(indv_list)):
                        if(indv_list[b] == info_ID):
                            if(merge_list[b + 4] == "."):
                                merge_list[b + 4] = info_list[a]
                            else:
                                merge_list[b + 4] = merge_list[b + 4] + "|" + info_list[a]
                        else:
                            continue
                
                # build merge string and output
                merged_string = "\t".join(merge_list)
                
                # SPLIT POSITIONS IF CNVS ARE MISMATCHED
                
                # determine if CNVs are mismatched                
                # declare variables outside inline loop
                inmerge_CNV = "ERROR"
                CNV_mismatch = False
                caller_mismatch = False
                
                # separate line into tab-delimited list
                line_list = merged_string.split("\t")
                
                # see how many versions of the inmerge_CNV are there
                CNV_list = []
                
                # iterate through line column-by-column
                for i in range(0, len(line_list)):
                    if("/" in line_list[i]): # upon finding a genotype filed
                        if(inmerge_CNV == "ERROR"): # replace placeholder with first CNV
                            inmerge_CNV = line_list[i].split("/")[0]
                            CNV_list.append(line_list[i].split("/")[0])
                        elif(inmerge_CNV == line_list[i].split("/")[0]): # if matching, proceed as normal
                            CNV_list.append(line_list[i].split("/")[0])
                            continue
                        else:
                            CNV_mismatch = True
                            CNV_list.append(line_list[i].split("/")[0])
                
                
                
                # if CNVs are mismatched, perform split
                if(CNV_mismatch == True):
                    # if it's only duplications, write out as normal. If it's not:
                    if("0" in CNV_list or "1" in CNV_list):        # are there deletions?    
                        if any(int(x) >= 2 for x in CNV_list): # then there are duplications and a split is necessary
                        
                            # call split method to return a given del_string and a given dup_string
                            del_string = (CNV_split(line_list, CNV_list, caller_mismatch))[0]
                            dup_string = (CNV_split(line_list, CNV_list, caller_mismatch))[1]
                            
                            # write out results of split method
                            file_out.write(del_string + "\n")
                            file_out.write(dup_string + "\n")
                            
                            # write out where callers differ on the same individual
                            if(caller_mismatch == True):
                                caller_mismatch_out.write(dup_string + "\n")
                            
                            # reset ongoing merge variables
                            merge_ongoing = False
                            caller_mismatch = False
                            qsnp_index = 0
                            penncnv_index = 0
                            final_outer_start = 0
                            final_outer_stop = 0
                            outer_start = 0
                            outer_stop = 0 
                            info_list = []      
                            
                                
                        # if it's just deletions, output as normal
                        else:
                            
                            del_list = merged_string.split("\t")
                            
                            for j in range(0, len(del_list)):
                                # check if callers differ in CNV               
                                if("|" in del_list[j]):
                                    del_indv_list = del_list[j].split("|")
                                    qsnp_CNV = int(del_indv_list[0][0])
                                    penncnv_CNV = int(del_indv_list[1][0])
                                    if(qsnp_CNV != penncnv_CNV):
                                        caller_mismatch = True
                                        
                                    # assign avg del
                                    if(caller_mismatch == True):
                                        for k in range(0, len(del_indv_list)):
                                            del_col_list = del_indv_list[k].split("/")
                                            avg_del = round(((qsnp_CNV + penncnv_CNV)/2),1)
                                            del_col_list[0] = str(avg_del)
                                            del_indv_list[k] = "/".join(del_col_list)
                                    del_list[j] = "|".join(del_indv_list)
                            
                            del_list[3] = "DELETION ONLY"
                            del_string = "\t".join(del_list)
                            file_out.write(del_string + "\n")
                        
                            # reset ongoing merge variables
                            merge_ongoing = False
                            caller_mismatch = False
                            qsnp_index = 0
                            penncnv_index = 0
                            final_outer_start = 0
                            final_outer_stop = 0
                            outer_start = 0
                            outer_stop = 0 
                            info_list = []    
                            
                            
                    # if it's just duplications, output as normal
                    else:
                        line_list = merged_string.split("\t")

                        # calculate average CNV for duplications
                        dup_CNVs = CNV_list.copy()
                        dup_CNVs[:] = [x for x in dup_CNVs if x != "0"]
                        dup_CNVs[:] = [x for x in dup_CNVs if x != "1"]
                        dup_CNVs = [ int(x) for x in dup_CNVs ]
                
                        dup_list = line_list

                        for j in range(0, len(dup_list)):
                            if("|" in dup_list[j]):
                                dup_indv_list = dup_list[j].split("|")
                                # check if callers differ in CNV
                                qsnp_CNV = int(dup_indv_list[0][0])
                                penncnv_CNV = int(dup_indv_list[1][0])
                                if(qsnp_CNV != penncnv_CNV):
                                    caller_mismatch = True
                                            
                                # write out where callers differ on the same individual
                                if(caller_mismatch == True):
                                    for k in range(0, len(dup_indv_list)):
                                        dup_col_list = dup_indv_list[k].split("/")
                                        avg_dup = round(((qsnp_CNV + penncnv_CNV)/2),1)
                                        dup_col_list[0] = str(avg_dup)
                                        dup_indv_list[k] = "/".join(dup_col_list)
                                dup_list[j] = "|".join(dup_indv_list)
                        
                        dup_list[3] = "DUPLICATION ONLY"
                        dup_string = "\t".join(dup_list)
                        
                        # write out where callers differ on the same individual
                        if(caller_mismatch == True):
                            caller_mismatch_out.write(dup_string + "\n")

                        file_out.write(dup_string + "\n")
                        
                    
                        
                        # reset ongoing merge variables
                        merge_ongoing = False
                        caller_mismatch = False
                        qsnp_index = 0
                        penncnv_index = 0
                        final_outer_start = 0
                        final_outer_stop = 0
                        outer_start = 0
                        outer_stop = 0 
                        info_list = []         
                
                # if there is no CNV mismatch, output as normal
                else:
                    file_out.write(merged_string + "\n")
                    
                    
                    # reset ongoing merge variables
                    merge_ongoing = False
                    caller_mismatch = False
                    qsnp_index = 0
                    penncnv_index = 0
                    final_outer_start = 0
                    final_outer_stop = 0
                    outer_start = 0
                    outer_stop = 0 
                    info_list = []         
                
   
            # if there is no ongoing merge, just output the current line
            else:
                unmerged_list = curr_line[0:3]
                unmerged_list.append(split_string)
                
                # add columns to unmerged_list
                for l in range(0, len(indv_list)):
                    unmerged_list.append(".")
                
                
                # add caller information to the corresponding individual column
                for b in range(0, len(indv_list)):
                    if(indv_list[b] == curr_ID):
                        unmerged_list[b + 4] = curr_geno
                    else:
                        continue
                
                # output as normal
                merged_string = "\t".join(unmerged_list)
                
                file_out.write(merged_string + "\n")
                
                # reset ongoing merge variables
                merge_ongoing = False
                qsnp_index = 0
                penncnv_index = 0
                final_outer_start = 0
                final_outer_stop = 0
                outer_start = 0
                outer_stop = 0

    
    # close input files
    file_in.close()
    indv_in.close()
    error_out.close()
    
    # close output file
    file_out.close()
    caller_mismatch_out.close()
    
    
    print("CNV_merge script complete.")
    
    # return count of merges
    return(merge_count)
    
# METHOD TO REMERGE WITHOUT CHANGING FORMAT
def CNV_remerge(input_file, output_label, cycle_count, error_file):
    """
    > cycles in the output from CNV_merge()
    > retains format but performs same merge operations
    > continues merging until no more merges can be performed
    """
    # keep track of how many cycles CNV_remerge has undergone
    print("CNV_remerge script is running cycle #" + str(cycle_count) + "...")
    
    # open files
    file_in = open(input_file, "r")
    file_out = open(output_label, "w")
    error_out = open(error_file, "w")
    
    # REPEAT OVERLAP CHECK
    # read in file as newline-separated list
    file_list = file_in.read().split("\n")
    
    # out-of-loop variables
    merge_ongoing = False # boolean that is used to enable multiple-line merges
    info_list = []
    final_outer_start = 0
    final_outer_stop = 0
    remerges_needed = 0
    
    # write out header
    file_out.write(file_list[0] + "\n")
    
    # run through the entire sorted_list, reading two lines at a time
    for i in range(1, len(file_list) - 1):
        
        # declare current line and next line as lists
        curr_line = file_list[i].split("\t")
        next_line = file_list[i + 1].split("\t")
        
        # skip any empty lines at the end
        if(len(next_line) <= 1 or len(curr_line) <= 1):
            continue
        
        # assign positions and determine length for current line
        curr_chrom = curr_line[0]
        curr_start = int(curr_line[1])
        curr_stop = int(curr_line[2])
        curr_length = curr_stop - curr_start
        curr_type = curr_line[3]
        
        # assign positions and determine length for next line
        next_chrom = next_line[0]
        next_start = int(next_line[1])
        next_stop = int(next_line[2])
        next_length = next_stop - next_start
        next_type = next_line[3]
        
        # skip sex chromosomes
        if curr_chrom == "X" or curr_chrom == "Y":
            continue
        
        # declare variable to hold information of whether it's a mismatch and what kind
        split_string = "NO MISMATCH"
        
        # CHECK IF THEY OVERLAP SUCH THAT THEY SHOULD BE MERGED
        overlap = False
        
        # must first be from the same chromosome AND be the same type of CNV!
        if(curr_chrom == next_chrom and curr_type == next_type):
            
            # determine overlap (inner) and outer positions
            if(curr_start >= next_start): # OVERLAP START
                overlap_start = curr_start # greater start is rightmost or inner start
                outer_start = next_start # lesser start is leftmost or outer start
            else:
                overlap_start = next_start # greater start is rightmost or inner start
                outer_start = curr_start # lesser start is leftmost or outer start
                
            if(curr_stop <= next_stop): # OVERLAP STOP
                overlap_stop = curr_stop # lesser stop is the leftmost or inner stop
                outer_stop = next_stop # greater stop is rightmost or outer stop
            else:
                overlap_stop = next_stop # lesser stop is the leftmost or inner stop
                outer_stop = curr_stop # greater stop is the rightmost or outer stop
            
            # adjust final_outer variables
            if(final_outer_start == 0 or final_outer_start >= outer_start):
                final_outer_start = outer_start
            else:
                final_outer_start = final_outer_start
            if(final_outer_stop == 0 or final_outer_stop <= outer_stop):
                final_outer_stop = outer_stop
            else:
                final_outer_stop = final_outer_stop
            
            # calculate overlap length
            overlap_length = overlap_stop - overlap_start
            
            # perform overlap 70% calculation
            curr_70 = curr_length * .7
            next_70 = next_length * .7
            if(overlap_length >=  curr_70 and overlap_length >= next_70):
                overlap = True
            else:
                overlap = False
        # if chromosomes do not match, no overlap is possible
        else:
            overlap = False
            
            
        # diagnostic step to troubleshoot merge of 152069090
        if(curr_stop >= 151840000 and curr_stop <= 152070000):
            error_out.write(str(curr_chrom) + "\t" + str(curr_start)  + "\t" +  str(curr_stop) + ", the overlap condition is " + str(overlap) + "\n")
            error_out.write("\t".join(next_line) + "\n")
            error_out.write("curr_start = " + str(curr_start) + "\n")
            error_out.write("curr_stop = " + str(curr_stop) + "\n")
            error_out.write("curr_length = " + str(curr_length) + "\n")
            error_out.write("curr_70 = " + str(curr_70) + "\n")
            error_out.write("next_start = " + str(next_start) + "\n")
            error_out.write("next_stop = " + str(next_stop) + "\n")
            error_out.write("next_length = " + str(next_length) + "\n")
            error_out.write("next_70 = " + str(next_70) + "\n")
            error_out.write("overlap_start = " + str(overlap_start) + "\n")
            error_out.write("overlap_stop = " + str(overlap_stop) + "\n")
            error_out.write("overlap_length = " + str(overlap_length) + "\n")
            # calculate overlap percentages
            if(curr_length >= overlap_length):
                overlap_with_curr = overlap_length/curr_length
            else:
                overlap_with_curr = curr_length/overlap_length
            if(next_length >= overlap_length):
                overlap_with_next = overlap_length/next_length
            else:
                overlap_with_next = next_length/overlap_length
            error_out.write("overlap_with_curr = " + str(overlap_with_curr) + "\n")
            error_out.write("overlap_with_next = " + str(overlap_with_next) + "\n\n")

        # IF OVERLAP IS TRUE, PERFORM MERGE
        if(overlap == True):
            remerges_needed += 1
            
            # perform merge
            merge_ongoing = True # initiate merge process
            
            # declare merge list for output
            merge_list = [curr_chrom, str(final_outer_start), str(final_outer_stop), curr_type]
            
            # sort through curr_line and next_line to see where they differ on individual data
            for A in range(4, len(curr_line)):
                if curr_line[A] != next_line[A]:
                    if curr_line[A] == ".": # overwrite null fields with data
                        merge_list.append(next_line[A])
                    elif next_line[A] == ".":
                        merge_list.append(curr_line[A])
                    else:
                        chimaera_string = curr_line[A] + "|" + next_line[A]
                        merge_list.append(chimaera_string)
                else:
                    merge_list.append(curr_line[A]) # copy over the matching field
                
        else: # IF OVERLAP IS FALSE, EITHER OUTPUT ONGOING MERGE OR CURRENT LINE
            
            # if there is an ongoing merge, output it
            if(merge_ongoing == True):   
                # build merge string and output
                merged_string = "\t".join(merge_list) + "\n"
                file_out.write(merged_string)

                # reset ongoing merge variables
                merge_ongoing = False
                final_outer_start = 0
                final_outer_stop = 0
                outer_start = 0
                outer_stop = 0 
                
            else:
                file_out.write(file_list[i] + "\n")
                
                # reset ongoing merge variables
                merge_ongoing = False
                final_outer_start = 0
                final_outer_stop = 0
                outer_start = 0

            
    print("Remerges needed: " + str(remerges_needed))
            
    
    # close files
    file_in.close()
    file_out.close()
    error_out.close()
    
    print("CNV_remerge has finished running cycle #" + str(cycle_count))

# METHOD TO PERFORM AN ADJACENCY-AGNOSTIC MERGE
def AJA_remerge(input_file, output_label, cycle_count, error_file):
    """
    > sorts though remerge output
    > takes in every position, then searches every OTHER position for a match
    > when a match is found for Variant A with Variant B:
        > writes out merged line in place of Variant A
        > marks Variant B to be skipped when printing non-overlapping lines
    """
    # keep track of how many cycles CNV_remerge has undergone
    print("AJA_remerge script is running cycle #" + str(cycle_count) + "...")
    
    # open files
    file_in = open(input_file, "r")
    file_out = open(output_label, "w")
    error_out = open(error_file, "w")
    
    # REPEAT OVERLAP CHECK
    # read in file as newline-separated list
    file_list = file_in.read().split("\n")
    
    # out-of-loop variables
    merge_ongoing = False # boolean that is used to enable multiple-line merges
    info_list = []
    final_outer_start = 0
    final_outer_stop = 0
    remerges_needed = 0
    previously_overlapped = [] # list to store positions that have already been merged and do not need to be listed again
    
    # write out header
    file_out.write(file_list[0] + "\n")
    
    # run through the entire sorted_list, reading two lines at a time
    for i in range(1, len(file_list) - 1):
        
        # declare current line and associated variables
        curr_line = file_list[i].split("\t")
        curr_chrom = curr_line[0]
        curr_start = int(curr_line[1])
        curr_stop = int(curr_line[2])
        curr_length = curr_stop - curr_start
        curr_type = curr_line[3]
        
        overlap = False
        
        # compare current line to every other line
        for j in range(1, len(file_list) -1):
            # skip the curr_line
            if file_list[j] == file_list[i]:
                continue
            else:
                # declare variant B's associated variables
                next_line = file_list[j].split("\t")
                next_chrom = next_line[0]
                next_start = int(next_line[1])
                next_stop = int(next_line[2])
                next_length = next_stop - next_start
                next_type = next_line[3]
                
                # declare positions to test for skips
                position_to_test = [str(curr_chrom), str(curr_start), str(curr_stop), curr_type]
                position_to_test_string = "\t".join(position_to_test)
                
                # skip any empty lines at the end
                if(len(next_line) <= 1 or len(curr_line) <= 1):
                    continue
            
                # skip sex chromosomes
                if curr_chrom == "X" or curr_chrom == "Y":
                    continue
                
                # skip positions in previously_overlapped
                if position_to_test_string in previously_overlapped:
                    continue
                
                # declare variable to hold information of whether it's a mismatch and what kind
                split_string = "NO MISMATCH"
                
                # CHECK IF THEY OVERLAP SUCH THAT THEY SHOULD BE MERGED
                # must first be from the same chromosome AND be the same type of CNV!
                if(curr_chrom == next_chrom and curr_type == next_type):
                    
                    # determine overlap (inner) and outer positions
                    if(curr_start >= next_start): # OVERLAP START
                        overlap_start = curr_start # greater start is rightmost or inner start
                        outer_start = next_start # lesser start is leftmost or outer start
                    else:
                        overlap_start = next_start # greater start is rightmost or inner start
                        outer_start = curr_start # lesser start is leftmost or outer start
                        
                    if(curr_stop <= next_stop): # OVERLAP STOP
                        overlap_stop = curr_stop # lesser stop is the leftmost or inner stop
                        outer_stop = next_stop # greater stop is rightmost or outer stop
                    else:
                        overlap_stop = next_stop # lesser stop is the leftmost or inner stop
                        outer_stop = curr_stop # greater stop is the rightmost or outer stop
                    
                    # adjust final_outer variables
                    if(final_outer_start == 0 or final_outer_start >= outer_start):
                        final_outer_start = outer_start
                    else:
                        final_outer_start = final_outer_start
                    if(final_outer_stop == 0 or final_outer_stop <= outer_stop):
                        final_outer_stop = outer_stop
                    else:
                        final_outer_stop = final_outer_stop
                    
                    # calculate overlap length
                    overlap_length = overlap_stop - overlap_start
                    
                    # perform overlap 70% calculation
                    curr_70 = curr_length * .7
                    next_70 = next_length * .7
                    if(overlap_length >=  curr_70 and overlap_length >= next_70):
                        overlap = True
                    else:
                        overlap = False
                # if chromosomes do not match, no overlap is possible
                else:
                    overlap = False
        
                # IF OVERLAP IS TRUE, PERFORM MERGE AND OUTPUT
                if(overlap == True):
                    remerges_needed += 1
                    
                    # declare merge list for output
                    merge_list = [str(curr_chrom), str(final_outer_start), str(final_outer_stop), curr_type]
                    
                    # sort through curr_line and next_line to see where they differ on individual data
                    for A in range(4, len(curr_line)):
                        if curr_line[A] != next_line[A]:
                            if curr_line[A] == ".": # overwrite null fields with data
                                merge_list.append(next_line[A])
                            elif next_line[A] == ".":
                                merge_list.append(curr_line[A])
                            else:
                                chimaera_string = curr_line[A] + "|" + next_line[A]
                                merge_list.append(chimaera_string)
                        else:
                            merge_list.append(curr_line[A]) # copy over the matching field
                    
                    # write out merged line in place of current line
                    file_out.write("\t".join(merge_list) + "\n")
                    position_to_skip = [str(next_chrom), str(next_start), str(next_stop), next_type]
                    position_to_skip_string = "\t".join(position_to_skip)
                    previously_overlapped.append(position_to_skip_string)
            
                    # diagnostic
                    error_out.write(str(curr_chrom) + "\t" + str(curr_start)  + "\t" +  str(curr_stop) + ", the overlap condition is " + str(overlap) + "\n")
                    error_out.write("\t".join(next_line) + "\n")
                    error_out.write("curr_start = " + str(curr_start) + "\n")
                    error_out.write("curr_stop = " + str(curr_stop) + "\n")
                    error_out.write("curr_length = " + str(curr_length) + "\n")
                    error_out.write("curr_70 = " + str(curr_70) + "\n")
                    error_out.write("next_start = " + str(next_start) + "\n")
                    error_out.write("next_stop = " + str(next_stop) + "\n")
                    error_out.write("next_length = " + str(next_length) + "\n")
                    error_out.write("next_70 = " + str(next_70) + "\n")
                    error_out.write("overlap_start = " + str(overlap_start) + "\n")
                    error_out.write("overlap_stop = " + str(overlap_stop) + "\n")
                    error_out.write("overlap_length = " + str(overlap_length) + "\n")
                    # calculate overlap percentages
                    if(curr_length >= overlap_length):
                        overlap_with_curr = overlap_length/curr_length
                    else:
                        overlap_with_curr = curr_length/overlap_length
                    if(next_length >= overlap_length):
                        overlap_with_next = overlap_length/next_length
                    else:
                        overlap_with_next = next_length/overlap_length
                    error_out.write("overlap_with_curr = " + str(overlap_with_curr) + "\n")
                    error_out.write("overlap_with_next = " + str(overlap_with_next) + "\n\n")
            
                    # reset merge variables
                    final_outer_start = 0
                    final_outer_stop = 0
                    outer_start = 0
                    outer_stop = 0
                    
                    # break the j loop
                    break
                else:
                    # reset merge variables
                    final_outer_start = 0
                    final_outer_stop = 0
                    outer_start = 0
                    outer_stop = 0
                    continue
                
            # reset merge variables
            final_outer_start = 0
            final_outer_stop = 0
            outer_start = 0
            outer_stop = 0
        
        # if it has progressed through every other line and overlap is still false, then no known overlaps
        else:
            if position_to_test_string in previously_overlapped:
                continue
            
                # reset merge variables
                final_outer_start = 0
                final_outer_stop = 0
                outer_start = 0
                outer_stop = 0
                
            else:
                file_out.write(file_list[i] + "\n")
                
                # reset merge variables
                final_outer_start = 0
                final_outer_stop = 0
                outer_start = 0
                outer_stop = 0
            
    print("Remerges needed: " + str(remerges_needed))
            
    
    # close files
    file_in.close()
    file_out.close()
    error_out.close()
    
    print("AJA_remerge has finished running cycle #" + str(cycle_count))

CNV_merge("../Data/merge_outputs/axiom_sorted.tsv", "_merge_intermediary_1")
CNV_merge("../Data/merge_outputs/jan_sorted.tsv", "_merge_intermediary_1")
CNV_merge("../Data/merge_outputs/may_sorted.tsv", "_merge_intermediary_1")

# TO DO: Fix merging error
CNV_merge("../Data/merge_outputs/all_sorted.tsv", "_merge_intermediary_1")

# remerge consecutive positions
CNV_remerge("../Data/merge_outputs/all_merge_intermediary_1.tsv", "../Data/merge_outputs/all_merge_intermediary_2.tsv", 1, "error.txt") # returns 16 remerges
CNV_remerge("../Data/merge_outputs/all_merge_intermediary_2.tsv", "../Data/merge_outputs/all_merge_intermediary_3.tsv", 2, "error.txt") # returns 0 remerges

# adjacency-agnostic merge
AJA_remerge("../Data/merge_outputs/all_merge_intermediary_2.tsv", "../Data/merge_outputs/all_merge_intermediary_3.tsv", 1, "error.txt") # returns 120 remerges
AJA_remerge("../Data/merge_outputs/all_merge_intermediary_3.tsv", "../Data/merge_outputs/all_merge_intermediary_4.tsv", 2, "error.txt") # returns 29 remerges
AJA_remerge("../Data/merge_outputs/all_merge_intermediary_4.tsv", "../Data/merge_outputs/all_merge_intermediary_5.tsv", 3, "error.txt") # returns 8 remerges
AJA_remerge("../Data/merge_outputs/all_merge_intermediary_5.tsv", "../Data/merge_outputs/all_merge_intermediary_6.tsv", 4, "error.txt") # returns 2 remerges
AJA_remerge("../Data/merge_outputs/all_merge_intermediary_6.tsv", "../Data/merge_outputs/all_merge_intermediary_7.tsv", 5, "error.txt") # returns 1 remerges
AJA_remerge("../Data/merge_outputs/all_merge_intermediary_7.tsv", "../Data/merge_outputs/all_merge_final.tsv", 6, "error8.txt") # returns 0 remerges
