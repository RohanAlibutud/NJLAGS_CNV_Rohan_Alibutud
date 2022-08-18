"""
GMNP_dominant.py
> Takes in an annotated CNV tab-delimited file
> Calculates internal frequency of each variant and output it to a separate file
> Takes in Geminesque output
> Uses positions from Geminesque to search through the AnnotSV file
> Removes benign genotypes
> Checks against frequency file and removes common genotypes
> Unlike standard GMN_Prioritizer.py, does not require the "de novo" tag
"""

# INTERNAL FREQUENCY CALCULATOR METHOD
def find_frequencies(input_file, filename):
    
    # open files
    anno_in = open(input_file, "r")
    
    # declare dictionary to hold frequency data
    freq_dict = {}
    
    # split input file into lines and declare header
    anno_list = anno_in.read().split("\n")
    anno_list.pop(-1)
    header_list = anno_list[0].split("\t")
    
    # link individual IDs to column number
    indv_dict = {}
    for z in range(14, 537):
        indv_dict[str(z)] = header_list[z]
    
    # iterate through each variant
    for a in range(1, len(anno_list)):
        indv_list = anno_list[a].split("\t")
        
        # list for storing variant genotypes
        var_list = []
        # list for families represented by variants in individual
        fam_list = []
        
        # for each variant, find which individuals have alternative genotypes
        for b in range(14, 537):
            if indv_list[b] != "0/0:0:0:0:0:0:null:0":
                var_list.append(indv_list[b])
                
                # isolate individual ID to determine family relation
                indv_ID = indv_dict[str(b)]
                fam_ID = indv_ID[0:4]
                if fam_ID not in fam_list:
                    fam_list.append(fam_ID)
                else:
                    continue
                      
        # add frequency to dictionary
        var_pos = indv_list[0].split("_")
        pos_string = "_".join(var_pos[0:3])
        
        # change addition depending on whether it's present in more than one family
        if(len(fam_list) <= 1):
            freq_dict[pos_string] = str(len(var_list)) + "/619"
        else:
            freq_dict[pos_string] = str(len(var_list)) + "/619 and present in multiple families"
                
    # write out frequencies as dictionary to be opened later
    with open(filename + "frequencies.txt","w") as data: 
        data.write(str(freq_dict))
    
    # close files
    anno_in.close()
    
    # return the frequency dictionary
    return freq_dict


# VARIANT BASED PROCESSOR METHOD
def var_process(gmn_input, annotsv_input, output_filename):
    """
    var_process()
    > imports an annotated .tsv and Geminesque output
    > sorts through each variant and outputs relevant columns
    > incorporates data from Geminesque output
    """
    print("Running GMN_Prioritizer.py var_process() method...")

    # import files
    gmn_in = open(gmn_input, "r")
    annotsv_in = open(annotsv_input, "r")
    var_out = open(output_filename, "w")
    
    # read annotsv file into list and separate a header line 
    annotsv_list = annotsv_in.read().split("\n")
    header_list = annotsv_list[0].split("\t")
    
    # read Geminesque file into list
    gmn_list = gmn_in.read().split("\n")
    
    # write out output header line
    output_header_list = ["chrom", "start", "stop", "sv_type", "overlap_CDS_percent", "auto_dom_anno", "patho_gain_list[phen|hpo|source|coord]", "patho_loss_list[phen|hpo|source|coord]", "benign_gain_list[source|coord]", "benign_loss_list[source|coord]", "cohort_freq", "genes", "individuals"]
    output_header_string = "\t".join(output_header_list)
    var_out.write(output_header_string + "\n")
    
    # link individual IDs to column number
    indv_dict = {}
    for z in range(14, 537):
        indv_dict[str(z)] = header_list[z]
    
    # iterate through AnnotSV file
    for a in range(1, len(annotsv_list)-1):
        
        # retain relevant columns
        var_list = annotsv_list[a].split("\t")
        var_svid = var_list[0]
        sv_chrom = int(var_list[1])
        sv_start = int(var_list[2])
        sv_stop = int(var_list[3])
        sv_length = int(var_list[4])
        sv_type = var_list[5]
        gene_id = var_list[539]
        auto_dom_anno = ""
        anno_mode = var_list[538]
        # account for events wherein the overlap % column is empty:
        if(var_list[544] == ""):
            overlap_CDS_percent = "N/A"
        else:
            overlap_CDS_percent = int(var_list[545])
            
        pathogenic = "False"
        benign = "False"
        
        # GATHER INFORMATION FOR PATHOGENIC AND BENIGN ANNOTATION        
        # declare list for pathogenic gain columns
        p_gain_list = var_list[556:560]
        for x in range(0, len(p_gain_list)):
            if(p_gain_list[x] == ""):
                p_gain_list[x] = "."
        p_gain_string = "|".join(p_gain_list)
                
        # declare list for pathogenic loss columns
        p_loss_list = var_list[560:564]
        for x in range(0, len(p_loss_list)):
            if(p_loss_list[x] == ""):
                p_loss_list[x] = "."
        p_loss_string = "|".join(p_loss_list)
        
        # declare list for pathogenic ins columns
        p_ins_list = var_list[564:568]
        for x in range(0, len(p_ins_list)):
            if(p_ins_list[x] == ""):
                p_ins_list[x] = "."
        p_ins_string = "|".join(p_ins_list)
                
        
        # declare list for benign gain columns
        b_gain_list = var_list[570:572]
        for x in range(0, len(b_gain_list)):
            if(b_gain_list[x] == ""):
                b_gain_list[x] = "."
        b_gain_string = "|".join(b_gain_list)
                
        # declare list for benign loss columns
        b_loss_list = var_list[572:574]
        for x in range(0, len(b_loss_list)):
            if(b_loss_list[x] == ""):
                b_loss_list[x] = "."
        b_loss_string = "|".join(b_loss_list)
        
        # declare list for benign loss columns
        b_ins_list = var_list[574:576]
        for x in range(0, len(b_ins_list)):
            if(b_ins_list[x] == ""):
                b_ins_list[x] = "."
        b_ins_string = "|".join(b_ins_list)
        
        # compile list of individuals the variant is present in 
        indvs_present = []
        
        # calculate cohort frequency
        indv_list = []
        fam_list = []
        for c in range(14, 537):
            if var_list[c] != "0/0:0:0:0:0:0:null:0":
                indv_list.append(var_list[c])
                # isolate individual ID to determine family relation
                indv_ID = indv_dict[str(c)]
                indvs_present.append(indv_ID) # keep track of list of individuals the variant is present in
                fam_ID = indv_ID[0:4]
                if fam_ID not in fam_list:
                    fam_list.append(fam_ID)
                else:
                    continue
                
        # change addition depending on whether it's present in more than one family
        if(len(fam_list) <= 1):
            cohort_freq = str(len(indv_list)) + "/619"
        else:
            cohort_freq = str(len(indv_list)) + "/619 and present in multiple families"
            
        # add family info from the Geminesque file using SV_id as index
        for b in range(1, len(gmn_list)-2):
            gene_list = gmn_list[b].split("\t")
            gene_svid = gene_list[4]
            
            # compile list of individuals into a single column
            indvs_string = "|".join(indvs_present)
    
            # retain if SVID matches AND there is auto_dom annotation AND annotation mode is full
            auto_dom_field = gene_list[8]
            if(var_svid in gene_svid and auto_dom_field != "" and anno_mode == "full"):
                # assign Geminesque data
                auto_dom_anno = "family" + auto_dom_field.split("family")[1]
                # write out whole list
                output_list = [str(sv_chrom), str(sv_start), str(sv_stop), sv_type, str(overlap_CDS_percent), auto_dom_anno, p_gain_string, p_loss_string, b_gain_string, b_loss_string, cohort_freq, gene_id, indvs_string]
                output_string = "\t".join(output_list)
                var_out.write(output_string + "\n")
                break
            else: 
                continue 


    # close files
    gmn_in.close()
    annotsv_in.close()
    var_out.close()
    
    print("var_process() method has completed running.")

# CHECK OVERLAP SUBMETHOD
def check_overlap(sv_chrom, sv_start, sv_stop, test_chrom, test_start, test_stop):
    final_outer_start = 0
    final_outer_stop = 0
    test_length = test_stop - test_start
    sv_length = sv_stop - sv_start
    
    # must first be from the same chromosome
    if(sv_chrom == test_chrom):
        # determine overlap (inner) and outer positions
        if(test_start >= sv_start): # OVERLAP START
            overlap_start = test_start # greater start is rightmost or inner start
            outer_start = sv_start # lesser start is leftmost or outer start
        else:
            overlap_start = sv_start # greater start is rightmost or inner start
            outer_start = test_start # lesser start is leftmost or outer start
        if(test_stop <= sv_stop): # OVERLAP STOP
            overlap_stop = test_stop # lesser stop is the leftmost or inner stop
            outer_stop = sv_stop # greater stop is rightmost or outer stop
        else:
            overlap_stop = sv_stop # lesser stop is the leftmost or inner stop
            outer_stop = test_stop # greater stop is the rightmost or outer stop
        
        # adjust final_outer variables
        if(final_outer_start == 0 or final_outer_start >= sv_start):
            final_start = sv_start
        else:
            final_start = final_start
        if(final_outer_stop == 0 or final_outer_stop <= outer_stop):
            final_outer_stop = outer_stop
        else:
            final_outer_stop = final_outer_stop
        
        # calculate overlap length
        overlap_length = overlap_stop - overlap_start
        outer_length = outer_stop - outer_start
    
        # perform overlap 70% calculation
        test_70 = round(test_length * .7)
        sv_70 = round(sv_length * .7)
        if(overlap_length >= test_70):
            test_overlap = True
        else:
            test_overlap = False
    else:
        test_overlap = False

    return(test_overlap)

# CHROMOSOME VALUE CONVERTER SUBMETHOD
def chrom_value(input_list):
    """
    chrom_value
    > read in a chrom from the first column of a .vcf line
    > translate it into a numerical value for the sake of sorting
    """    
    
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
    
# VARIANT SORTER METHOD
def var_sort(var_input, var_output):
    """
    var_sort()
    > imports GMN_variants output and then sorts them
    > simple bubble sort
    """
    print("Running GMN_Prioritizer.py var_sort() method...")

    # import files
    var_file = open(var_input, "r")
    
    # declare export file
    var_out = open(var_output, "w")
    
    # PERFORM THE SORTING (bubble sort)
    unsorted_list = var_file.read().split("\n")
    
    # set holding list
    sorted_holding_list = []    
    
    # set sorting boolean
    swapped = True
    while swapped == True:
        swapped = False
        for i in range(1, len(unsorted_list) - 1):
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
                
    # append to sorted holding list
    for i in range(len(unsorted_list)):
         # skip any empty lines
        if(len(unsorted_list[i]) <= 1):
            continue
        
        # write to file
        sorted_holding_list.append(unsorted_list[i])
    
    # write to file
    sorted_string = "\n".join(sorted_holding_list)
    var_out.write(sorted_string)
        
    # close files
    var_file.close()
    var_out.close()

    print("var_sort() method has completed running.")

# VARIANT PRIORITIZER METHOD
def var_prioritize(input_file, output_file, disp_file, strictness):
    """
    var_prioritize()
    > imports sorted GMN variants and prioritizes them
    > sorts through each variant, filters/ranks
    > outputs passing variants
    """
    
    print("Running GMN_Prioritizer.py var_prioritize() method...")
    
    # open files
    file_in = open(input_file, "r")
    file_out = open(output_file + "_" + strictness, "w")
    disp_in = open(disp_file, "r")
    
    # build list of variants
    var_list = file_in.read().split("\n")
    
    # build list of dispensable genes
    disp_genes = disp_in.read().split("\n")
    
    # output header
    header_string = var_list[0] + "\t" + "dispensable"
    file_out.write(header_string + "\n")
    
    # iterate through list of variants and remove ones with matching benign annotation and no corresponding pathogenic
    for variant in range(1, len(var_list)):
        line_list = var_list[variant].split("\t")
        
        # define columns
        chrom = line_list[0]
        start = line_list[1]
        stop = line_list[2]
        sv_type = line_list[3]
        overlap_CDS = line_list[4]
        auto_dom = line_list[7]
        p_gain = line_list[6]
        p_loss = line_list[7]
        b_gain = line_list[8]
        b_loss = line_list[9]
        cohort_freq = line_list[10]
        genes = line_list[11]
        individuals = line_list[12]
        
        # FILTER ON PRESENCE OF DISPENSABLE GENES
        gene_list = genes.split(";") 
        
        disp_bool = False
        for gene in range(0, len(gene_list)):
            if(gene_list[gene] in disp_genes):
                disp_bool = True
            else:
                disp_bool = False
                break
        line_list.append(str(disp_bool))
        
        var_list[variant] = "\t".join(line_list)
        
        # FILTER ON PATHOGENIC/BENIGN ANNOTATION
        if strictness == "strict":
            # filter on presence of benign annotation
            if b_gain != ".|." and b_loss != ".|.":
                continue
            else:
                file_out.write(var_list[variant] + "\n")
        else:     
            # filter on presence of benign annotation AND absence of pathogenic annotation
            if b_gain != ".|." and b_loss != ".|." and p_gain == ".|.|.|." and p_loss == ".|.|.|.":
                continue
            else:
                file_out.write(var_list[variant] + "\n")
                
        
    
    # close files
    file_in.close()
    
    print("var_prioritize() method has completed running.")

# ADD OVERLAP_CDS METHOD
def add_overlap_CDS(var_input, CNV_input, var_output):
    """
    add_overlap_cds()
    > imports sorted and prioritized variants
    > replaces blank overlap_cds column with the list of corresponding values from the split columns
    > exports modified variant table

    """
    
    print("Running GMN_Prioritizer.py add_overlap_cds() method...")
    
    # open files
    file_in = open(var_input, "r")
    file_out = open(var_output, "w")
    CNV_in = open(CNV_input, "r")
    
    # build list of variants
    var_list = file_in.read().split("\n")
    # build list of variants from annotation file
    CNV_list = CNV_in.read().split("\n")
    
    # output header
    header_string = var_list[0]
    file_out.write(header_string + "\n")
    
    # iterate through list of variants and identify full lines as keys, and append split lines to them as values
    pos_list = []
    search_start = 1
    for variant in range(1, len(var_list)-1):
        line_list = var_list[variant].split("\t")
        
        # identify full variants by position
        chrom = line_list[0]
        start = line_list[1]
        stop = line_list[2]
        var_pos = str(chrom) + "_" + str(start) + "_" + str(stop)
        
        # append to pos list
        pos_list.append(var_pos)
        
        # declare holding list for overlap_CDS values
        oCDS_list = []
        
        # use var_pos to identify the matching split lines in the CNV2_anno file
        merging_splits = False
        for a in range(search_start, len(CNV_list)-1):
            anno_list = CNV_list[a].split("\t")
                
            # identify position currently being examined
            anno_pos = "_".join((anno_list[0].split("_")[0:3]))
            
            # check annotation mode
            anno_mode = anno_list[538]
            
            # compare var_pos to anno_pos
            if(var_pos == anno_pos and anno_mode == "full"):
                merging_splits = True         
                continue
            elif(var_pos == anno_pos and anno_mode == "split"):
                oCDS_list.append(anno_list[546])
            elif(merging_splits == True):
                # export oCDS list for this variant
                line_list[4] = ":".join(oCDS_list)
                merging_splits = False
                break
        
        # export variant
        file_out.write("\t".join(line_list) + "\n")
                
       
        
    
    print("add_overlap_CDS() method has completed running.")

# ADD GTEx EXPRESSION DATA METHOD
def add_gtex(var_input, var_output, gtex_input):
    """
    add_gtex()
    > imports sorted and prioritized variants
    > searches all the genes listed from each variant in the gtex file and returns tissue expression data
    > exports modified variant table

    """
    
    print("Running GMN_Prioritizer.py add_gtex() method...")
    
    # open files
    file_in = open(var_input, "r")
    file_out = open(var_output, "w")
    gtex_in = open(gtex_input, "r")
    
    # declare lists
    file_list = file_in.read().split("\n")
    gtex_list = gtex_in.read().split("\n")
    
    # declare header
    header_list = file_list[0].split("\t")
    header_list.append("GTEx [Amygdala | Anterior cingulate cortex (BA24) | Caudate (basal ganglia) | Cerebellar Hemisphere | Cerebellum | Cortex | Frontal Cortex (BA9) | Hippocampus | Hypothalamus | Nucleus accumbens (basal ganglia) | Putamen (basal ganglia) | Spinal cord (cervical c-1) | Substantia nigra]")
    header_string = "\t".join(header_list)
    file_out.write(header_string + "\n")
    
    # sort through each of the variants and extract all genes
    for variant in range(1, len(file_list)-1):
        var_list = file_list[variant].split("\t")
        genes = var_list[11].split(";")
        
        # declare GTEx output list
        gtex_out_list = []
        
        # sort through all the genes and use them as keys to search through the GTEx data
        for gene in range(0, len(genes)):     
            
            
            # search through GTEx file
            for a in range(1, len(gtex_list)):
                
                # isolate columns for brain expression
                brain_cols = gtex_list[a].split("\t")[9:22]
                brain_cols = [round(float(num), 2) for num in brain_cols]
                brain_cols = [str(num) for num in brain_cols]
                
                # if the gene is found in the gtex_list, print out its expression data for all tissues
                if(genes[gene] in gtex_list[a]):
                    gtex_string = genes[gene] + " : " + ";".join(brain_cols)
                    gtex_out_list.append(gtex_string)
            
        # append gtex out list to variant output
        gtex_out_string = "|".join(gtex_out_list)
        var_list.append(gtex_out_string)
        var_out_string = "\t".join(var_list)
        file_out.write(var_out_string + "\n")
                
    # close files
    file_in.close()
    file_out.close()
    gtex_in.close()
    
    print("add_gtex() method has completed running.")

# ADD EXPRESSION DATA FROM ALL THREE SOURCES METHOD
def add_expression(var_input, var_output, tpm1, tpm2, tpm3):
    """
    add_expression()
    > takes in variant-based .tsv input
    > references three tissue expression databases
    > requires every variant to have at least one gene with tpm > 3 for a given brain tissue
    > exports passing variants
    """
    print("Running add_expression() method...")

    import pandas as pd
    import pickle
    
    # PROCESS TPM1
    df_tpm1 = pd.read_csv(tpm1, sep = "\t")
    columns_GTEx = list(df_tpm1.columns[1:-1]) # select GTEx columns from dataframe
    columns_GTEx_brain = [e for e in columns_GTEx if "Brain" in e] # retain only the ones from brain 
    df_tpm1['tpm1_max'] = df_tpm1[columns_GTEx_brain].max(axis=1) # declare new column that collects maximum tpm value from all brain columns
    # reduce df_tpm1 to only the relevant columns
    df_tpm1.rename(columns={'Description': 'GENE'}, inplace=True)
    df_tpm1 = df_tpm1[['GENE', 'tpm1_max']]
    
    
    # PROCESS TPM2
    df_tpm2 = pd.read_csv(tpm2, sep = "\t")
    # treat all columns as necessary for now
    columns_brainspan = list(df_tpm2.columns[2:]) 
    df_tpm2['tpm2_max'] = df_tpm2[columns_brainspan].max(axis=1)

    
    # PROCESS TPM3
    df_tpm3 = pickle.load(open(tpm3,'rb'))#165 tissues 
    df_tpm3['tpm3_max'] = df_tpm3[df_tpm3.columns[2:]].max(axis=1) 
    df_tpm3 = df_tpm3[['Gene_Name', 'tpm3_max'] + [e for e in df_tpm3.columns if "cerebral_cortex" in e]]
    # The second column is what we care about: TPM3
    df_tpm3.rename(columns={'Gene_Name': 'GENE'}, inplace = True)
    df_tpm3 = df_tpm3[['GENE', 'tpm3_max']]
    

    # declare separate dataframe to hold all values from the various tpm files
    df_tpm_1and2 = df_tpm1.merge(df_tpm2, how = 'left', left_on = "GENE", right_on ='geneSymbol')
    df_tpm_all = df_tpm_1and2.merge(df_tpm3, how = 'left', left_on = "GENE", right_on ='GENE') 
    
    
    
    # COMPARE ALL VARIANTS TO DF_TPM_ALL
    # open variant file and export file
    var_in = open(var_input, "r")
    var_list = var_in.read().split("\n")
    var_out = open(var_output, "w")
    
    # write out header list
    header_list = var_list[0].split("\t")
    header_list.append("Expression")
    var_out.write("\t".join(header_list) + "\n")
    
    # cast df_tpm_all to dict
    tpm_dict = df_tpm_all.set_index("GENE").T.to_dict("list")
    
    # iterate through variant file line by line
    for a in range(1, len(var_list)-1):
        line_list = var_list[a].split("\t")
        gene_list = line_list[11].split(";")
        
        # boolean to determine if any tpm exceeds notability threshold (>3 or >5)
        is_expr = False
        
        # iterate through every gene in variant, and if 0 genes have enough expression, drop it
        for b in range(0, len(gene_list)):
            gene_id = gene_list[b]
            
            if gene_id not in tpm_dict:
                continue
            else:
                gene_tpm = tpm_dict[gene_id]
                # iterate through all the tpm information for a given gene
                for c in range(3, len(gene_tpm)):
                    tpm_value = round(float(gene_tpm[c]), 2)
                    if tpm_value >= 3:
                        is_expr = True
                        break
                    else:
                        continue
        
    
        # after having read through all the genes in the variant, decide to write out or not:
        """
        if is_expr == True:
            var_out.write(var_list[a] + "\n")
        else:
            continue
        """
        var_out.write(var_list[a] + "\t" + str(is_expr) + "\n")
          
    print("add_expression() method has completed.")
            

# COMMANDS
phenotype = "LI_star"
filepath = "/lab01/Projects/Rohan_Projects/CNV_Project/2022/Total_Pathway/Data/auto_dom_analysis/"
var_process(filepath + phenotype + "/" + phenotype + "_doms.tsv", filepath + "CNV2_anno_corrected.tsv",filepath + phenotype + "/GMNP_intermediate1.tsv")
var_sort(filepath + phenotype + "/GMNP_intermediate1.tsv", filepath + phenotype + "/GMNP_intermediate2.tsv")
var_prioritize(filepath + phenotype + "/GMNP_intermediate2.tsv", filepath + phenotype + "/GMNP_intermediate3.tsv", filepath + "dispensable_genes_and_muc.txt", "lenient")
add_overlap_CDS(filepath + phenotype + "/GMNP_intermediate3.tsv_lenient", filepath + "CNV2_anno_corrected.tsv", filepath + phenotype + "/GMNP_intermediate4.tsv")
add_gtex(filepath + phenotype + "/GMNP_intermediate4.tsv", filepath + phenotype + "/prioritized_variants.tsv", filepath + "GTEx_2019_12_12.txt")


add_expression(filepath + phenotype + "/GMNP_intermediate4.tsv", filepath + phenotype + "/GMNP_intermediate5.tsv", filepath + "tpm1.txt", filepath + "tpm2.txt", filepath + "tpm3.txt")
