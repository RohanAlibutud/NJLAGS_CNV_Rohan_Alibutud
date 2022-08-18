# add in OVerlapped_CDS_percent column from AnnotSV
# define function to hold process
def add_Overlapped():
    print("Running Overlapped() method")
    
    # open relevant files
    GMN_in = open("../Data/GMN_out_2022_01_14.tsv", "r")
    AnnotSV_in = open("../Data/CNV2_anno_corrected.tsv", "r")
    GMN_out = open("../Data/GMN_out_2022_01_14_final.tsv", "w")

    # read out geneIDs from the GMN_out file
    GMN_list = GMN_in.read().split("\n")
    anno_list = AnnotSV_in.read().split("\n")
    
    # insert new information to header
    header_list = GMN_list[0].split("\t")
    header_list.insert(3,"Overlapped_CDS_percent")
    header_string = "\t".join(header_list)
    
    # write header out to file    
    GMN_out.write(header_string + "\n")
    
    # iterate through GMN_list and add AnnotSV info
    for i in range(1, len(GMN_list)):
        line_list = GMN_list[i].split("\t")
        var_gene_ID = line_list[0]        
    
        # use geneIDs to search through AnnotSV 
        for j in range(1, len(anno_list)):
            # only search through lines in the "split" annotation mode
            if("split" not in anno_list[j]):
                continue
            else:
                row_list = anno_list[j].split("\t")                
                anno_gene_ID = row_list[539] # geneID
                if(var_gene_ID == anno_gene_ID):
                    anno_OCDS = row_list[546]
                    break
        
        line_list.insert(3, anno_OCDS)
        output_string = "\t".join(line_list)
        GMN_out.write(output_string + "\n")
        
        
        
    
    # close files
    GMN_in.close()
    AnnotSV_in.close()
    GMN_out.close()

print("Overlap method called")
add_Overlapped()
print("Overlap method completed")

