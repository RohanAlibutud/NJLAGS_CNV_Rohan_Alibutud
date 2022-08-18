# all the following scripts require os
import os
import statistics


"""
CNV_combine
> take in .csv files from cnv calls folder
> combine without merging or sorting
> output a combined .tsv file
"""    
def CNV_combine(input_folder):
    print("CNV_combine script running...")
    
    # MERGE AXIOM BATCH
    # open axiom files
    axiom_qsnp_file = open(input_folder + "axiom_run3_combined_20200306_lbf_v2_bpfix_info.csv", "r")
    axiom_penncnv_file = open(input_folder + "axiom_run7_20200302_gc_all_v2_bpfix_info.csv", "r")    
    axiom_unsorted_file = open("../Data/merge_outputs/axiom_unsorted.tsv", "w")
    
    # declare merged list for both axiom callers
    axiom_merge_list = []
    
    # define column positions for axiom_qsnp_file
    for line in axiom_qsnp_file:
        line_list = line.split(",")
        for i,j in enumerate(line_list):
            if(line_list[i] == "partid"):
                partid_loc = i
            elif(line_list[i] == "chromosome" or line_list[i] == "chrom"):
                old_chrom_loc = i
                chrom_loc = i
            elif(line_list[i] == "start position (bp)" or line_list[i] == "start" or line_list[i] == "start_pos"):
                old_start_loc = i
                start_loc = i
            elif(line_list[i] == "end position (bp)" or line_list[i] == "stop" or line_list[i] == "stop_pos"):
                old_stop_loc = i
                stop_loc = i
            elif(line_list[i] == "copy number" or line_list[i] == "cnv"):
                copy_loc = i
            elif(line_list[i] == "no. probes"  or line_list[i] == "numsnp"):
                probes_loc = i
            elif(line_list[i] == "cnv qc"):
                cnv_qc_loc = i
            else:
                continue
        break
    
    # import variants from axiom qsnp file
    axiom_qsnp_list = [] # axiom qsnp holding list
    for line in axiom_qsnp_file:
        
        # extract positions from the header of that file
        line_list = line.split(",")
        header_len = len(line_list)
            
        # assign chrom, start, stop, partid
        chrom = line_list[chrom_loc]
        start = line_list[start_loc]
        stop = line_list[stop_loc]
        partid = line_list[partid_loc]
        
        # assign qnsp information
        old_chrom = line_list[old_chrom_loc]
        old_start = line_list[old_start_loc]
        old_stop = line_list[old_stop_loc]
        copy_num = line_list[copy_loc]
        cnv_qc = line_list[cnv_qc_loc] 
        probes = line_list[probes_loc]
        qsnp_list = ["qsnp", old_chrom, old_start, old_stop, copy_num, cnv_qc, probes]
        qsnp_string = ":".join(qsnp_list)
            
        # define output list
        if(chrom == "23"):
            chrom = "X"
        if(chrom == "24"):
            chrom = "Y"
        output_list = [chrom, start, stop, partid, qsnp_string]
        output_string = "\t".join(output_list)
        # skip any empty lines
        if(len(output_string) <= 1):
            continue
        
        axiom_qsnp_list.append(output_string)
    
    # turn axiom_list into a string
    axiom_qsnp_string = "\n".join(axiom_qsnp_list)
        
    # write to file
    axiom_unsorted_file.write(axiom_qsnp_string)
    
    # separate qsnp from penncnv
    axiom_unsorted_file.write("\n")
        
    # define column positions for axiom_penncnv_file
    for line in axiom_penncnv_file:
        line_list = line.split(",")
        for i,j in enumerate(line_list):
            if(line_list[i] == "partid"):
                partid_loc = i
            elif(line_list[i] == "chromosome" or line_list[i] == "chrom"):
                old_chrom_loc = i
                chrom_loc = i
            elif(line_list[i] == "start position (bp)" or line_list[i] == "start" or line_list[i] == "start_pos"):
                old_start_loc = i
                start_loc = i
            elif(line_list[i] == "end position (bp)" or line_list[i] == "stop" or line_list[i] == "stop_pos"):
                old_stop_loc = i
                stop_loc = i
            elif(line_list[i] == "copy number" or line_list[i] == "cnv"):
                copy_loc = i
            elif(line_list[i] == "no. probes" or line_list[i] == "numsnp"):
                probes_loc = i
            elif(line_list[i] == "cnv qc"):
                cnv_qc_loc = i
            else:
                continue
        break
    
    # import variants from axiom penncnv file    
    axiom_penncnv_list = [] # axiom penncnv holding list
    for line in axiom_penncnv_file:
        
        # extract positions from the header of that file
        line_list = line.split(",")
        header_len = len(line_list)
            
        # assign chrom, start, stop, partid
        chrom = line_list[chrom_loc]
        start = line_list[start_loc]
        stop = line_list[stop_loc]
        partid = line_list[partid_loc]
        
        # assign qnsp information
        old_chrom = line_list[old_chrom_loc]
        old_start = line_list[old_start_loc]
        old_stop = line_list[old_stop_loc]
        copy_num = line_list[copy_loc]
        cnv_qc = line_list[cnv_qc_loc] 
        probes = line_list[probes_loc].strip()
        penncnv_list = ["penncnv", old_chrom, old_start, old_stop, copy_num, cnv_qc, probes]
        penncnv_string = ":".join(penncnv_list)
            
        # define output list
        if(chrom == "X"):
            continue
        if(chrom == "Y"):
            continue
        output_list = [chrom, start, stop, partid, penncnv_string]
        output_string = "\t".join(output_list)
        
        # skip any empty lines
        if(len(output_string) <= 1):
            continue
        
        # append to holding list
        axiom_penncnv_list.append(output_string)
    
    # join axiom_penncnv_list to string
    axiom_penncnv_string = "\n".join(axiom_penncnv_list)
        
    # write to file
    axiom_unsorted_file.write(axiom_penncnv_string)
        
    # close intermediate file
    axiom_unsorted_file.close()
    
    # close axiom files
    axiom_qsnp_file.close()
    axiom_penncnv_file.close()
    
    
    
    # MERGE JAN 2015 ILLUMINA BATCH
    # open jan 2015 files
    jan_qsnp_file = open(input_folder + "illm_jan2015_run2_combined_20200306_lbf_v2_bpfix_info.csv", "r")
    jan_penncnv_file = open(input_folder + "psych_run4_20200302_gc_fixed_all_v2_bpfix_info.csv", "r")    
    jan_unsorted_file = open("../Data/merge_outputs/jan_unsorted.tsv", "w")
    
    # declare merged list for both axiom callers
    jan_merge_list = []
    
    # define column positions for 
    for line in jan_qsnp_file:
        line_list = line.split(",")
        for i,j in enumerate(line_list):
            if(line_list[i] == "new_stop"):
                stop_loc = i
            elif(line_list[i] == "partid"):
                partid_loc = i
            elif(line_list[i] == "chromosome" or line_list[i] == "chrom"):
                old_chrom_loc = i
                chrom_loc = i
            elif(line_list[i] == "start position (bp)" or line_list[i] == "start" or line_list[i] == "start_pos"):
                old_start_loc = i
                start_loc = i
            elif(line_list[i] == "end position (bp)" or line_list[i] == "stop" or line_list[i] == "stop_pos"):
                old_stop_loc = i
                stop_loc = i
            elif(line_list[i] == "copy number" or line_list[i] == "cnv"):
                copy_loc = i
            elif(line_list[i] == "no. probes" or line_list[i] == "numsnp"):
                probes_loc = i
            elif(line_list[i] == "cnv qc"):
                cnv_qc_loc = i
            else:
                continue
        break
    
    # import variants from jan qsnp file  
    jan_qsnp_list = [] # jan qsnp holding list
    for line in jan_qsnp_file:
        
        # extract positions from the header of that file
        line_list = line.split(",")
        header_len = len(line_list)
            
        # assign chrom, start, stop, partid
        chrom = line_list[chrom_loc]
        start = line_list[start_loc]
        stop = line_list[stop_loc]
        partid = line_list[partid_loc]
        
        # assign qnsp information
        old_chrom = line_list[old_chrom_loc]
        old_start = line_list[old_start_loc]
        old_stop = line_list[old_stop_loc]
        copy_num = line_list[copy_loc]
        cnv_qc = line_list[cnv_qc_loc] 
        probes = line_list[probes_loc]
        qsnp_list = ["qsnp", old_chrom, old_start, old_stop, copy_num, cnv_qc, probes]
        qsnp_string = ":".join(qsnp_list)
            
        # define output list
        if(chrom == "23"):
            chrom = "X"
        if(chrom == "24"):
            chrom = "Y"
        output_list = [chrom, start, stop, partid, qsnp_string]
        output_string = "\t".join(output_list)
        # skip any empty lines
        if(len(output_string) <= 1):
            continue
        
        # jan qsnp list
        jan_qsnp_list.append(output_string)
    
    # join list into outputabble string
    jan_qsnp_string = "\n".join(jan_qsnp_list)
    
    # write to file
    jan_unsorted_file.write(jan_qsnp_string)
        
        
    # define column positions for jan_penncnv_file
    for line in jan_penncnv_file:
        line_list = line.split(",")
        for i,j in enumerate(line_list):
            if(line_list[i] == "partid"):
                partid_loc = i
            elif(line_list[i] == "chromosome" or line_list[i] == "chrom"):
                old_chrom_loc = i
                chrom_loc = i
            elif(line_list[i] == "start position (bp)" or line_list[i] == "start" or line_list[i] == "start_pos"):
                old_start_loc = i
                start_loc = i
            elif(line_list[i] == "end position (bp)" or line_list[i] == "stop" or line_list[i] == "stop_pos"):
                old_stop_loc = i
                stop_loc = i
            elif(line_list[i] == "copy number" or line_list[i] == "cnv"):
                copy_loc = i
            elif(line_list[i] == "no. probes" or line_list[i] == "numsnp"):
                probes_loc = i
            elif(line_list[i] == "cnv qc"):
                cnv_qc_loc = i
            else:
                continue
        break
    
    # separate qsnp from penncnv
    jan_unsorted_file.write("\n")
    
    # import variants from axiom penncnv file    
    jan_penncnv_list = [] # jan penncnv holding list
    for line in jan_penncnv_file:
        
        # extract positions from the header of that file
        line_list = line.split(",")
        header_len = len(line_list)
            
        # assign chrom, start, stop, partid
        chrom = line_list[chrom_loc]
        start = line_list[start_loc]
        stop = line_list[stop_loc]
        partid = line_list[partid_loc]
        
        # assign qnsp information
        old_chrom = line_list[old_chrom_loc]
        old_start = line_list[old_start_loc]
        old_stop = line_list[old_stop_loc]
        copy_num = line_list[copy_loc]
        cnv_qc = line_list[cnv_qc_loc] 
        probes = line_list[probes_loc].strip()
        penncnv_list = ["penncnv", old_chrom, old_start, old_stop, copy_num, cnv_qc, probes]
        penncnv_string = ":".join(penncnv_list)
            
        # define output list
        if(chrom == "X"):
            continue
        if(chrom == "Y"):
            continue
        output_list = [chrom, start, stop, partid, penncnv_string]
        output_string = "\t".join(output_list)
        
         # skip any empty lines
        if(len(output_string) <= 1):
            continue
        
        # append output string to holding list
        jan_penncnv_list.append(output_string)
        
    # join holding list into a string
    jan_penncnv_string = "\n".join(jan_penncnv_list)
     
    # write to file
    jan_unsorted_file.write(jan_penncnv_string)
        
    
    # close intermediate file
    jan_unsorted_file.close()
    
    
    # close axiom files
    jan_qsnp_file.close()
    jan_penncnv_file.close()
    
    
    
    # MERGE MAY 2017 ILLUMINA BATCH
    # open may 2017 files
    may_qsnp_file = open(input_folder + "illm_may2017_run3_combined_20200306_lbf_v2_bpfix_info.csv", "r")
    may_penncnv_file = open(input_folder + "psych_may2017_run3_20200302_gc_all_v2_bpfix_info.csv", "r")    
    may_unsorted_file = open("../Data/merge_outputs/may_unsorted.tsv", "w")
    
    # declare merged list for both may callers
    may_merge_list = []
    
    # define column positions for may_qsnp_file
    for line in may_qsnp_file:
        line_list = line.split(",")
        for i,j in enumerate(line_list):
            if(line_list[i] == "partid"):
                partid_loc = i
            elif(line_list[i] == "chromosome" or line_list[i] == "chrom"):
                old_chrom_loc = i
                chrom_loc = i
            elif(line_list[i] == "start position (bp)" or line_list[i] == "start" or line_list[i] == "start_pos"):
                old_start_loc = i
                start_loc = i
            elif(line_list[i] == "end position (bp)" or line_list[i] == "stop" or line_list[i] == "stop_pos"):
                old_stop_loc = i
                stop_loc = i
            elif(line_list[i] == "copy number" or line_list[i] == "cnv"):
                copy_loc = i
            elif(line_list[i] == "no. probes" or line_list[i] == "numsnp"):
                probes_loc = i
            elif(line_list[i] == "cnv qc"):
                cnv_qc_loc = i
            else:
                continue
        break
    
    # import variants from may qsnp file    
    may_qsnp_list = []
    for line in may_qsnp_file:
        
        # extract positions from the header of that file
        line_list = line.split(",")
        header_len = len(line_list)
            
        # assign chrom, start, stop, partid
        chrom = line_list[chrom_loc]
        start = line_list[start_loc]
        stop = line_list[stop_loc]
        partid = line_list[partid_loc]
        
        # assign qnsp information
        old_chrom = line_list[old_chrom_loc]
        old_start = line_list[old_start_loc]
        old_stop = line_list[old_stop_loc]
        copy_num = line_list[copy_loc]
        cnv_qc = line_list[cnv_qc_loc] 
        probes = line_list[probes_loc]
        qsnp_list = ["qsnp", old_chrom, old_start, old_stop, copy_num, cnv_qc, probes]
        qsnp_string = ":".join(qsnp_list)
            
        # define output list
        if(chrom == "23"):
            chrom = "X"
        if(chrom == "24"):
            chrom = "Y"
        output_list = [chrom, start, stop, partid, qsnp_string]
        output_string = "\t".join(output_list)
        # skip any empty lines
        if(len(output_string) <= 1):
            continue
        
        # append output string to holding list
        may_qsnp_list.append(output_string)
    
    # join holding list into output string
    may_qsnp_string = "\n".join(may_qsnp_list)
        
    # write to file
    may_unsorted_file.write(may_qsnp_string)
    
    # separate qsnp from penncnv
    may_unsorted_file.write("\n")
        
    # define column positions for may_penncnv_file
    for line in may_penncnv_file:
        line_list = line.split(",")
        for i,j in enumerate(line_list):
            if(line_list[i] == "partid"):
                partid_loc = i
            elif(line_list[i] == "chromosome" or line_list[i] == "chrom"):
                old_chrom_loc = i
                chrom_loc = i
            elif(line_list[i] == "start position (bp)" or line_list[i] == "start" or line_list[i] == "start_pos"):
                old_start_loc = i
                start_loc = i
            elif(line_list[i] == "end position (bp)" or line_list[i] == "stop" or line_list[i] == "stop_pos"):
                old_stop_loc = i
                stop_loc = i
            elif(line_list[i] == "copy number" or line_list[i] == "cnv"):
                copy_loc = i
            elif(line_list[i] == "no. probes" or line_list[i] == "numsnp"):
                probes_loc = i
            elif(line_list[i] == "cnv qc"):
                cnv_qc_loc = i
            else:
                continue
        break
    
    # import variants from may penncnv file
    may_penncnv_list = []
    for line in may_penncnv_file:
        
        # extract positions from the header of that file
        line_list = line.split(",")
        header_len = len(line_list)
            
        # assign chrom, start, stop, partid
        chrom = line_list[chrom_loc]
        start = line_list[start_loc]
        stop = line_list[stop_loc]
        partid = line_list[partid_loc]
        
        # assign qnsp information
        old_chrom = line_list[old_chrom_loc]
        old_start = line_list[old_start_loc]
        old_stop = line_list[old_stop_loc]
        copy_num = line_list[copy_loc]
        cnv_qc = line_list[cnv_qc_loc] 
        probes = line_list[probes_loc].strip()
        penncnv_list = ["penncnv", old_chrom, old_start, old_stop, copy_num, cnv_qc, probes]
        penncnv_string = ":".join(penncnv_list)
            
        # define output list
        if(chrom == "X"):
            continue
        if(chrom == "Y"):
            continue
        output_list = [chrom, start, stop, partid, penncnv_string]
        output_string = "\t".join(output_list)
        
        # skip any empty lines
        if(len(output_string) <= 1):
            continue
        
        # append output_string to holding list
        may_penncnv_list.append(output_string)
        
    # join holding list into output string
    may_penncnv_string = "\n".join(may_penncnv_list)
        
    # write to file
    may_unsorted_file.write(may_penncnv_string)
        
    
    # close intermediate file
    may_unsorted_file.close()
    
    
    # close may files
    may_qsnp_file.close()
    may_penncnv_file.close()
    
    print("Combining complete.")

"""
chrom_value
> read in a chrom from the first column of a .vcf line
> translate it into a numerical value for the sake of sorting
"""    
# CHROMOSOME VALUE CONVERTER
def chrom_value(input_list):
    if(input_list == ""):
        return(0)
    
    input_list = input_list.split("\t")
    
    chrom_val = input_list[0]
    if(chrom_val == "X"):
        chrom_val = 24
    elif(chrom_val == "Y"):
        chrom_val = 25
    else:
        chrom_val = int(chrom_val)
    return(chrom_val)

"""
CNV_sort
> take in a unsorted .tsv file
> use a bubble sort to order calls by chromosome and positions
> output a sorted .tsv file
"""    
def CNV_sort(input_file):
    print("CNV_sort script running...")
    
    # parse input filename
    if("axiom" in input_file):
        filehandle = "axiom"
    elif("jan" in input_file):
        filehandle = "jan"
    elif("may" in input_file):
        filehandle = "may"
    elif("all" in input_file):
        filehandle = "all"
    
    # open input and output files
    unsorted_file = open(input_file, "r")
    sorted_file = open(projects_folder + "/Data/merge_outputs/" + filehandle + "_sorted.tsv", "w")
    
    # PERFORM THE SORTING (bubble sort)
    unsorted_list = unsorted_file.read().split("\n")
    
    # set holding list
    sorted_holding_list = []    
    
    # set sorting boolean
    swapped = True
    while swapped == True:
        swapped = False
        for i in range(len(unsorted_list) - 1):
            # skip empty lines
            if(unsorted_list[i] == "\n"):
                continue
            
            # break up each line in unsorted_list into columns
            line1_list = unsorted_list[i].split("\t")
            line2_list = unsorted_list[i + 1].split("\t")
            
            # sort by chromosome
            if(chrom_value(unsorted_list[i]) > chrom_value(unsorted_list[i + 1])):
                # swap the elements
                unsorted_list[i], unsorted_list[i + 1] = unsorted_list[i + 1], unsorted_list[i]
                swapped = True
            # sort by start position
            elif((chrom_value(unsorted_list[i]) == chrom_value(unsorted_list[i + 1])) and (int(line1_list[1]) > int(line2_list[1]))):
                # swap the elements
                unsorted_list[i], unsorted_list[i + 1] = unsorted_list[i + 1], unsorted_list[i]
                swapped = True
            # sort by stop position
            elif((chrom_value(unsorted_list[i]) == chrom_value(unsorted_list[i + 1])) and (int(line1_list[1]) == int(line2_list[1])) and (int(line1_list[2]) > int(line2_list[2]))):
                # swap the elements
                unsorted_list[i], unsorted_list[i + 1] = unsorted_list[i + 1], unsorted_list[i]
                swapped = True
            # sort by individual
            elif((chrom_value(unsorted_list[i]) == chrom_value(unsorted_list[i + 1])) and (int(line1_list[1]) == int(line2_list[1])) and (int(line1_list[2]) == int(line2_list[2])) and (int(line1_list[3]) > int(line2_list[3]))):    
                # swap the elements
                unsorted_list[i], unsorted_list[i + 1] = unsorted_list[i + 1], unsorted_list[i]
                swapped = True
                
    # append to sorted holding list
    for i in range(len(unsorted_list)):
         # skip any empty lines
        if(len(unsorted_list[i]) <= 1):
            continue
        
        # write to file
        sorted_holding_list.append(unsorted_list[i])
    
    # write to file
    sorted_string = "\n".join(sorted_holding_list)
    sorted_file.write(sorted_string)
        
    # close files
    unsorted_file.close()
    sorted_file.close()
    
    print("Sorting complete.")
    
    
"""
CNV_merge
> take in a sorted .tsv file
> merge adjacent calls together
> output a merged .tsv file
"""    
def CNV_merge(input_file):
    print("CNV_merge script running...")
    
    # parse input filename
    if("axiom" in input_file):
        filehandle = "axiom"
    elif("jan" in input_file):
        filehandle = "jan"
    elif("may" in input_file):
        filehandle = "may"
    elif("test" in input_file):
        filehandle = "test"
    
    # open input and output files
    sorted_file = open(projects_folder + "/Data/merge_outputs/" + filehandle + "_sorted.tsv", "r")
    merged_file = open(projects_folder + "/Data/merge_outputs/" + filehandle + "_merged.tsv", "w")
    
    # import information in the sorted file as a list
    sorted_list = sorted_file.read().split("\n")
    info_list = [] # used for storing an increasing amount of info between loops
    qsnp_index = 0
    penncnv_index = 0
    
    # print out header
    header_list = ["chrom", "start", "stop", "partid", "caller_info[caller, chrom, start, stop, copy_number, qc, probes]\n"]
    header_string = "\t".join(header_list)
    merged_file.write(header_string)
    
    # out-of-loop variables
    merge_ongoing = False # boolean that is used to enable multiple-line merges
    
    # run through the entire sorted_list, reading two lines at a time
    for i in range(0, len(sorted_list) - 1):
        
        # declare current line and next line as lists
        curr_line = sorted_list[i].split("\t")
        next_line = sorted_list[i + 1].split("\t")
        
        # skip any empty lines at the end
        if(len(next_line) <= 1 or len(curr_line) <= 1):
            continue
        
        # assign positions and determine length for current line
        curr_chrom = curr_line[0]
        curr_start = int(curr_line[1])
        curr_stop = int(curr_line[2])
        curr_ID = curr_line[3]
        curr_info = curr_line[4]
        curr_length = curr_stop - curr_start
        
        # assign positions and determine length for next line
        next_chrom = next_line[0]
        next_start = int(next_line[1])
        next_stop = int(next_line[2])
        next_ID = next_line[3]
        next_info = next_line[4]
        next_length = next_stop - next_start
        
        # CHECK IF THEY OVERLAP SUCH THAT THEY SHOULD BE MERGED
        overlap = False        
        
        # must first be from the same individual and chromosome
        if(curr_ID == next_ID and curr_chrom == next_chrom):
            
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
            
            # calculate overlap length
            overlap_length = overlap_stop - overlap_start
            
            # perform overlap 70% calculation
            curr_70 = curr_length * .7
            next_70 = next_length * .7
            if(overlap_length >=  curr_70 and overlap_length >= next_70):
                overlap = True
            else:
                overlap = False
        
        # if chromosome and individual do not match, no overlap is possible
        else:
            overlap = False
            
        # if overlap is true, perform merge
        if(overlap == True):
            
            # PERFORM MERGE
            merge_ongoing = True # initiate merge process
            
            # declare merge list for output
            merge_list = [curr_chrom, str(outer_start), str(outer_stop), curr_ID]
            
            # append caller information
            curr_info_list = curr_info.split(":")
            if("qsnp" in curr_info_list[0]):
                qsnp_index += 1
                curr_info_list[0] = curr_info_list[0] + "_" + str(qsnp_index)
            elif("penncnv" in curr_info_list[0]):
                penncnv_index += 1
                curr_info_list[0] = curr_info_list[0] + "_" + str(penncnv_index)
            curr_info = ":".join(curr_info_list)
            info_list.append(curr_info)
            
        else:
            # if there is an ongoing merge, output it first
            if(merge_ongoing == True):                
                # append caller information for last line in merge
                curr_info_list = curr_info.split(":")
                if("qsnp" in curr_info_list[0]):
                    qsnp_index += 1
                    curr_info_list[0] = curr_info_list[0] + "_" + str(qsnp_index)
                elif("penncnv" in curr_info_list[0]):
                    penncnv_index += 1
                    curr_info_list[0] = curr_info_list[0] + "_" + str(penncnv_index)
                curr_info = ":".join(curr_info_list)
                
                # append caller info
                info_list.append(curr_info)
                info_string = "|".join(info_list)
                
                # build merge string and output
                merge_list.append(info_string)
                merge_string = "\t".join(merge_list)
                merged_file.write(merge_string + "\n")
                
                # reset ongoing merge variables
                merge_ongoing = False
                qsnp_index = 0
                penncnv_index = 0
                info_list = []
                
                
            # if there is no ongoing merge, just output the current line
            else:
                merged_file.write("\t".join(curr_line) + "\n")
        
    # account for last line if unmerged
    if(merge_ongoing == True):
        pass
    else:
        merged_file.write(sorted_list[-1])
    
            
    # close input and output files
    sorted_file.close()
    merged_file.close()
    print("Merging complete.")
    
    
"""
CNV_report
> takes in merged .tsv files
> takes counts and summary statistics
> outputs a summary statistics file
"""
def CNV_report(input_file):
    print("CNV_report script running...")
    
    # open input and output files
    file_in = open(input_file, "r")
    file_out = open(input_file.strip(".tsv") + "_report.txt", "w")
    data_out = open(input_file.strip(".tsv") + "_sum_stats.txt", "w")
    
    # declare count variables
    both_callers_count = 0
    qsnp_call_count = 0
    penncnv_call_count = 0
    identical_overlaps = 0
    nonidentical_overlaps = 0
    
    # declare list variables
    indv_list = []
    call_lengths_list = []
    
    # count how many of each caller show up
    for line in file_in:
        if("chrom" in line):
            continue
        
        # divide into columns for reading
        line_list = line.split("\t")
        indv_ID = line_list[3]
        
        # divide info section into columns
        if("|" in line):
            info_list = line_list[4].split("|")
            if("qsnp" in info_list[0]):
                qsnp_list = info_list[0].split(":")
                penncnv_list = info_list[1].split(":")
            else:
                qsnp_list = info_list[1].split(":")
                penncnv_list = info_list[0].split(":")
        else:
            info_list = line_list[4].split(":")
        
        # reset booleans
        both_bool = False
        qsnp_bool = False
        penncnv_bool = False    
        identical_start = False
        identical_stop = False
        
        # count callers
        if("qsnp" in line and "penncnv" in line):
            both_callers_count += 1
            both_bool = True
        elif("qsnp" in line):
            qsnp_call_count += 1
            qsnp_bool = True
        elif("penncnv" in line):
            penncnv_call_count += 1
            penncnv_bool = True
        else:
            print("ERROR")
            print(line)
        
        # count individuals
        if(indv_ID not in indv_list):
            indv_list.append(indv_ID)
        else:
            pass    
        
        # compare type of overlap
        if(both_bool == True):
            # are the positions identical or slightly different?
            if(qsnp_list[2] == penncnv_list[2]):
                identical_start = True
            if(qsnp_list[3] == penncnv_list[3]):
                identical_stop = True
            if(identical_start == True and identical_stop == True):
                identical_overlaps += 1
            else:
                nonidentical_overlaps += 1
                #print(line)
        
        # collate lengths of each call
        call_length = int(line_list[2]) - int(line_list[1])
        call_lengths_list.append(call_length)
    
    # calculate average length of a call
    average_call_length = statistics.median(call_lengths_list)
    # find longest call
    call_lengths_list.sort()
    longest_call = call_lengths_list[-1]
    # find shortest call
    shortest_call = call_lengths_list[0]
    
    # export the lengths of calls to data_out
    for i,j in enumerate(call_lengths_list, start = 0):
        data_out.write(str(call_lengths_list[i]) + "\n")
    
    # write outputs
    file_out.write("Calls present in both: " + str(both_callers_count) + "\n")
    file_out.write("Calls present in qsnp only: " + str(qsnp_call_count) + "\n")
    file_out.write("Calls present in penncnv only: " + str(penncnv_call_count) + "\n\n")
    file_out.write("Unique individuals present in callset: " + str(len(indv_list)) + "\n\n")
    file_out.write("Identical overlaps: " + str(identical_overlaps) + "\n")
    file_out.write("Nonidentical overlaps: " + str(nonidentical_overlaps) + "\n\n")
    file_out.write("Average length of a call: " + str("{:.2f}".format(average_call_length)) + "\n")
    file_out.write("Length of the longest call: " + str(longest_call) + "\n")
    file_out.write("Length of the shortest call: " + str(shortest_call) + "\n")
    
    # close input and output files
    file_in.close()
    file_out.close()
    data_out.close()
    
    print("Report finished.")
    

"""
CNV_report2
> takes in merged .tsv files
> takes counts and summary statistics
> outputs a summary statistics file
"""
def CNV_report2(input_file):
    print("CNV_report script running...")
    
    # open input and output files
    file_in = open(input_file, "r")
    file_out = open(input_file.strip(".tsv") + "_report.txt", "w")
    data_out = open(input_file.strip(".tsv") + "_sum_stats.txt", "w")
    
    # declare count variables
    both_callers_count = 0
    qsnp_call_count = 0
    penncnv_call_count = 0
    identical_overlaps = 0
    nonidentical_overlaps = 0
    
    # declare list variables
    indv_list = []
    call_lengths_list = []
    
    # count how many of each caller show up
    for line in file_in:
        if("chrom" in line):
            continue
        
        # divide into columns for reading
        line_list = line.split("\t")
        indv_ID = line_list[3]
        
        # divide info section into columns
        if("|" in line):
            info_list = line_list[4].split("|")
            if("qsnp" in info_list[0]):
                qsnp_list = info_list[0].split(":")
                penncnv_list = info_list[1].split(":")
            else:
                qsnp_list = info_list[1].split(":")
                penncnv_list = info_list[0].split(":")
        else:
            info_list = line_list[4].split(":")
        
        # reset booleans
        both_bool = False
        qsnp_bool = False
        penncnv_bool = False    
        identical_start = False
        identical_stop = False
        
        # count callers
        if("qsnp" in line and "penncnv" in line):
            both_callers_count += 1
            both_bool = True
        elif("qsnp" in line):
            qsnp_call_count += 1
            qsnp_bool = True
        elif("penncnv" in line):
            penncnv_call_count += 1
            penncnv_bool = True
        else:
            print("ERROR")
            print(line)
        
        # count individuals
        if(indv_ID not in indv_list):
            indv_list.append(indv_ID)
        else:
            pass    
        
        # compare type of overlap
        if(both_bool == True):
            # are the positions identical or slightly different?
            if(qsnp_list[2] == penncnv_list[2]):
                identical_start = True
            if(qsnp_list[3] == penncnv_list[3]):
                identical_stop = True
            if(identical_start == True and identical_stop == True):
                identical_overlaps += 1
            else:
                nonidentical_overlaps += 1
                #print(line)
        
        # collate lengths of each call
        call_length = int(line_list[2]) - int(line_list[1])
        call_lengths_list.append(call_length)
    
    # calculate average length of a call
    average_call_length = statistics.median(call_lengths_list)
    # find longest call
    call_lengths_list.sort()
    longest_call = call_lengths_list[-1]
    # find shortest call
    shortest_call = call_lengths_list[0]
    
    # export the lengths of calls to data_out
    for i,j in enumerate(call_lengths_list, start = 0):
        data_out.write(str(call_lengths_list[i]) + "\n")
    
    # write outputs
    file_out.write("Calls present in both: " + str(both_callers_count) + "\n")
    file_out.write("Calls present in qsnp only: " + str(qsnp_call_count) + "\n")
    file_out.write("Calls present in penncnv only: " + str(penncnv_call_count) + "\n\n")
    file_out.write("Unique individuals present in callset: " + str(len(indv_list)) + "\n\n")
    file_out.write("Identical overlaps: " + str(identical_overlaps) + "\n")
    file_out.write("Nonidentical overlaps: " + str(nonidentical_overlaps) + "\n\n")
    file_out.write("Average length of a call: " + str("{:.2f}".format(average_call_length)) + "\n")
    file_out.write("Length of the longest call: " + str(longest_call) + "\n")
    file_out.write("Length of the shortest call: " + str(shortest_call) + "\n")
    
    # close input and output files
    file_in.close()
    file_out.close()
    data_out.close()
    
    print("Report finished.")


projects_folder = "/lab01/Projects/Rohan_Projects/CNV_Project/2022/Total_Pathway"

CNV_combine(projects_folder + "/Data/starting_call_files/")

#CNV_sort(projects_folder + "/Data/merge_outputs/axiom_unsorted.tsv")
#CNV_merge(projects_folder + "/Data/merge_outputs/axiom_sorted.tsv")
   
#CNV_sort(projects_folder + "/Data/merge_outputs/jan_unsorted.tsv")
#CNV_merge(projects_folder + "/Data/merge_outputs/jan_sorted.tsv")

#CNV_sort(projects_folder + "/Data/merge_outputs/may_unsorted.tsv")
#CNV_merge(projects_folder + "/Data/merge_outputs/may_sorted.tsv")

# merge to write and form all_sorted
axiom_in = open(projects_folder + "/Data/merge_outputs/axiom_unsorted.tsv", "r")
jan_in = open(projects_folder + "/Data/merge_outputs/jan_unsorted.tsv", "r")
may_in = open(projects_folder + "/Data/merge_outputs/may_unsorted.tsv", "r")
all_out = open(projects_folder + "/Data/merge_outputs/all_unsorted.tsv", "w")

axiom_string = axiom_in.read()
jan_string = jan_in.read()
may_string = may_in.read()

all_string = axiom_string + "\n" + jan_string + "\n" + may_string
all_out.write(all_string)

axiom_in.close()
jan_in.close()
may_in.close()
all_out.close()

CNV_sort(projects_folder + "/Data/merge_outputs/all_unsorted.tsv")

