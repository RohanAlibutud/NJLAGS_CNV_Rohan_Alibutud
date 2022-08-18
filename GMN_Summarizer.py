"""
GMN_Summarizer.py
> Read in the tab-delimited Geminesque output file
> Iterate through the data and collate summary statistics
> Output to a .txt file
"""

# SUMMARY STATS METHOD
def gmn_summarize(input_file):
    print("gmn_summarize() is running on " + input_file)
    
    # declare dictionary for counting stats
    sum_dict = {
        "Gene count": 0,
        "De novo": 0,
        "Autosomal recessive": 0,
        "Autosomal dominant": 0,
        "Multiple segregation types": 0,
        "Multiple variants of same segregation type": 0,
        "Variants in CDSaffected": 0}
    
    # collate individual stats
    mul_var_genes = []
    indv_list = []

    # open files
    gmn_in = open(input_file, "r")
    #input_file.strip(".tsv") + "_stat.tsv"
    sum_out = open("GMN_out_FINAL_stats.tsv", "w")
    
    # import and iterate through gmn_in
    gmn_list = gmn_in.read().split("\n")    
    
    for i in range(1, len(gmn_list)-2):
        gene_list = gmn_list[i].split("\t")
        sum_dict["Gene count"] += 1
        
        # count variant types
        if(gene_list[6] != ""):
            sum_dict["De novo"] += 1
        if(gene_list[7] != ""):
            sum_dict["Autosomal recessive"] += 1
        if(gene_list[8] != ""):
            sum_dict["Autosomal dominant"] += 1
            
        # determine if a gene has multiple segregation types
        if(gene_list[6] != "" and gene_list[7] != "" or gene_list[7] != "" and gene_list[8] != "" or gene_list[6] != "" and gene_list[8] != ""):
            sum_dict["Multiple segregation types"] += 1
        
        # determine if gene has multiple variants of same segregation type of either type
        mul_var = False
        if(gene_list[6].count("variant_pos") > 1):
            mul_var = True
        if(gene_list[7].count("variant_pos") > 1):
            mul_var = True      
        if(gene_list[8].count("variant_pos") > 1):
            mul_var = True      
        if(mul_var == True):
            sum_dict["Multiple variants of same segregation type"] += 1
            if(gene_list[6] != ""):
                mul_var_genes.append(gene_list[0] + " [De novo]: " + gene_list[6])
            if(gene_list[7] != ""):
                mul_var_genes.append(gene_list[0] + " [Auto rec]: " + gene_list[7])
            if(gene_list[8] != ""):
                mul_var_genes.append(gene_list[0] + "[Auto dom]: " + gene_list[8])
                
        # which affected individual had the most variants?    
        for j in range(0, len(gene_list)):
            if("variant_pos" in gene_list[j]):
                var_list = gene_list[j].split(",")
                for k in range(0, len(var_list)):
                    if("family" in var_list[k]):
                        fam_list = var_list[k].split(":")
                        indv_ID = fam_list[2].split("(")[0]
                        indv_list.append(indv_ID)
                        
        # check if the gene's variants overlapped the coding region
        if(gene_list[2] == "True"):
            sum_dict["Variants in CDSaffected"] += 1
            
    
    # build dictionary of individuals
    indv_list = list(dict.fromkeys(indv_list))
    indv_dict = {}
    for l in range(0,len(indv_list)):
        indv_dict[indv_list[l]] = 0
        
    # count which individuals had the most variants
    # most de novo variants
    for m in range(1, len(gmn_list)-2):
        gene_list = gmn_list[m].split("\t")
        if("variant_pos" in gene_list[6]):
            var_list = gene_list[6].split(",")
            for o in range(0, len(var_list)):
                if("family" in var_list[o]):
                    fam_list = var_list[o].split(":")
                    indv_ID = fam_list[2].split("(")[0]
                    for p in range(0, len(indv_list)):
                        if(indv_list[p] == indv_ID):
                            indv_dict[indv_list[p]] += 1
    de_novo_dict = indv_dict
    dn_list = sorted(de_novo_dict.items(), key=lambda x: x[1], reverse=True)
    
    # most auto rec variants
    # reset
    indv_dict = {}
    for l in range(0,len(indv_list)):
        indv_dict[indv_list[l]] = 0
    
    # count
    for m in range(1, len(gmn_list)-2):
        gene_list = gmn_list[m].split("\t")
        if("variant_pos" in gene_list[7]):
            var_list = gene_list[7].split(",")
            for o in range(0, len(var_list)):
                if("family" in var_list[o]):
                    fam_list = var_list[o].split(":")
                    indv_ID = fam_list[2].split("(")[0]
                    for p in range(0, len(indv_list)):
                        if(indv_list[p] == indv_ID):
                            indv_dict[indv_list[p]] += 1
    auto_rec_dict = indv_dict
    ar_list = sorted(auto_rec_dict.items(), key=lambda x: x[1], reverse = True)
    
    
    # write to output file
    # basic stats
    for key, value in sum_dict.items():
        sum_tup = (key, ": ", str(value),"\n")
        sum_string = "".join(sum_tup)
        sum_out.write(sum_string)
    sum_out.write("\n")
    
    # list of genes with Multiple variants of same segregation type
    sum_out.write("Genes with multiple variants of same segregation type: \n")
    sum_out.write("\n".join(mul_var_genes) + "\n\n")
    
    # list of individuals by gene count
    sum_out.write("List of individuals with most de novo variants in genes, minimum 1:\n")
    for variant in dn_list:
        if(variant[1] > 0):
            sum_out.write(variant[0] + ": " + str(variant[1]) + "\n")
    sum_out.write("\n")
            
    sum_out.write("List of individuals with most autosomal recessive variants in genes, minimum 1:\n")
    for variant in ar_list:
        if(variant[1] > 0):
            sum_out.write(variant[0] + ": " + str(variant[1]) + "\n")
    
    # close files
    gmn_in.close()
    sum_out.close()
    
    print("gmn_summarize() is finished running.")

# DERIVED STATS METHOD
def gmn_derive(input_file):
    print("gmn_derive() is running on " + input_file)
    
    # open files
    gmn_in = open(input_file, "r")
    #input_file.strip(".tsv") + "_stat.tsv"
    der_out = open("GMN_out_FINAL_stats.tsv", "a")
    
    # declare variables for derived stats
    CNV_lengths = []
    de_novo_CNV_lengths = []
    auto_rec_CNV_lengths = []
    overlap_list = []
    dn_overlap_list = []
    ar_overlap_list = []
    highly_interesting = []
    
    # import and iterate through gmn_in
    gmn_list = gmn_in.read().split("\n")    
    for i in range(1, len(gmn_list)-2):
        gene_list = gmn_list[i].split("\t")
        
        # calculate length of the variants
        pos_list = gene_list[4].split("_")
        chrom = pos_list[0]
        start = pos_list[1]
        stop = pos_list[2]
        length = int(stop) - int(start)
        CNV_lengths.append(length)
        # check for segregation patterns
        if(gene_list[6] != ""):
            de_novo_CNV_lengths.append(length)
        if(gene_list[7] != ""):
            auto_rec_CNV_lengths.append(length)
        
        # calculate average percentage overlap with CDS
        overlap = int(gene_list[3])
        overlap_list.append(overlap)
        # check for segregation patterns
        if(gene_list[6] != ""):
            dn_overlap_list.append(overlap)
        if(gene_list[7] != ""):
            ar_overlap_list.append(overlap)
            
        # qualify gene for ranking
        if(gene_list[6] != "" and overlap >= 100 and length >= 1000000): # gene is de novo and mostly in coding region and long
            highly_interesting.append([length, gmn_list[i]])
            
        
    # perform calculations
    mean_length = round(sum(CNV_lengths)/len(CNV_lengths))
    mean_dn_length = round(sum(de_novo_CNV_lengths)/len(de_novo_CNV_lengths))
    mean_ar_length = round(sum(auto_rec_CNV_lengths)/len(auto_rec_CNV_lengths))
    mean_overlap = round(sum(overlap_list)/len(overlap_list), 2)
    mean_dn_overlap = round(sum(dn_overlap_list)/len(dn_overlap_list), 2)
    mean_ar_overlap = round(sum(ar_overlap_list)/len(ar_overlap_list), 2)
    
    # write outputs
    der_out.write("\n")
    der_out.write("Mean length of all CNVs: " + str(mean_length) + "\n")
    der_out.write("Mean length of all de novo CNVs: " + str(mean_dn_length) + "\n")
    der_out.write("Mean length of all autosomal recessive CNVS: " + str(mean_ar_length) + "\n")
    der_out.write("Mean coding region overlap percentage of all CNVs: " + str(mean_overlap) + "\n")
    der_out.write("Mean coding region overlap percentage of all de novo CNVs: " + str(mean_dn_overlap) + "\n")
    der_out.write("Mean coding region overlap percentage of all autosomal recessive CNVs: " + str(mean_ar_overlap) + "\n")
    
    # write out highly interesting genes
    der_out.write("\n")
    der_out.write("Highly interesting genes: \n")
    hi_dict = {}
    for gene in enumerate(highly_interesting):
        hi_dict[gene[0]] = gene[1]
    hi_dict_sorted = sorted(hi_dict.items(), key=lambda x: x[1], reverse=True)
    for i in range(0, len(hi_dict_sorted)):
        der_out.write(str(hi_dict_sorted[i][1][0]) + "  |  " + hi_dict_sorted[i][1][1] + "\n")
    
        
    # close files
    gmn_in.close()
    der_out.close()
    
    print("gmn_derive() has finished running.")
    
# COMMANDS
gmn_summarize("./GMN_out_FINAL.tsv")
gmn_derive("./GMN_out_FINAL.tsv")