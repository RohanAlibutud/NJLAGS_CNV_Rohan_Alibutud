"""
Manipulator.py
> reads through a file as tab-delimited table
> edits each cell as necessary
> returns as an edited file
"""

# COMPARE .PED AND .VCF
def vcf_ped_compare(input_file):
    file_in = open(input_file, "r")
    
    vcf_list = []
    ped_list = []
    
    for line in file_in:
        line_list = line.split("\t")
        vcf_list.append(line_list[0])
        ped_list.append(line_list[1].strip("\n"))

    for i in range(0, len(vcf_list)):
        if(vcf_list[i] not in ped_list and vcf_list[i] != ""):
            print(vcf_list[i])
        
    file_in.close()
        

# MANIPULATE A .PED FILE
def ped_manipulate(input_file, output_file):
    # import file
    file_in = open(input_file, "r")
    file_out = open(output_file, "w")
    
    # iterate throughout the file
    for line in file_in:
        line_list = line.strip("\n").split("\t")
        
        # identify columns
        family_id = line_list[0]
        part_id = line_list[1]
        wgs_id = line_list[2]
        ru_id = line_list[3]
        paternal_id = line_list[4]
        maternal_id = line_list[5]
        sex = line_list[6]
        ASD = line_list[8]
        
        # skip family 2018
        if(family_id == "2018"):
            continue
        
        # write out to output file
        file_out.write(line)
    
    # close files
    file_in.close()
    file_out.close()
    

# MANIPULATE A .VCF FILE
def vcf_manipulate(input_file, output_file):
    # open files
    file_in = open(input_file, "r")
    file_out = open(output_file, "w")
    
    # declare list for positions to be removed
    to_remove = []
    
    # iterate through every line of the file
    for line in file_in:
        # skip past specification info
        if(line[1] == "#"):
            file_out.write(line)
            continue
        
        # identify index position of family to remove
        if(line[0] == "#"):
            # divide header into list by tab
            header_list = line.strip("\n").split("\t")
        
            # iterate through header_list to find desired column
            for a in range(0, len(header_list)):
                if("2018" in header_list[a]):
                    to_remove.append(a)
                if("2117001" in header_list[a]):
                    to_remove.append(a)
                if("2108" in header_list[a]):
                    to_remove.append(a)
                if("2078" in header_list[a]):
                    to_remove.append(a)
                if("2094" in header_list[a]):
                    to_remove.append(a)
                if("2111" in header_list[a]):
                    to_remove.append(a)
                if("2119003" in header_list[a]):
                    to_remove.append(a)
            
            # remove target column in header
            for b in range(0, len(to_remove)):
                header_list.pop(to_remove[b]-b)
            
            # write out header
            header_string = "\t".join(header_list)
            file_out.write(header_string + "\n")
        else:
            # divide line into list by tab
            line_list = line.strip("\n").split("\t")
            
            # remove target columns in body
            for c in range(0, len(to_remove)):
                line_list.pop(to_remove[c]-c)
            
            # write out edited line
            line_string = "\t".join(line_list)
            file_out.write(line_string + "\n")
    
    # close files
    file_in.close()
    file_out.close()
    
# MANIPULATE AN ANNOTATED .TSV FILE
def tsv_manipulate(input_file, output_file):
    
    # import
    file_in = open(input_file, "r")
    file_out = open(output_file, "w")
    
    # iterate through file line by line
    for line in file_in:
        line_list = line.split("\t")
        
        
        #iterate through line
        for a in range(9, len(line_list)):
            # replace all nulls with homozygous reference
            geno_list = line_list[a].split(":")
            if(geno_list[0] == "./."):
                geno_list[0] = "0/0"
            # and reassign genotypes according to indel counts
            else:
                # 0.5/. should be called as 0/1 for heterozygous deletion
                if(geno_list[0] == "0.5/."):
                    geno_list[0] = "0/1"
                # 1.0/. should be called as 0/1 for heterozygous deletion
                elif(geno_list[0] == "1.0/."):
                    geno_list[0] = "0/1"
                # 3.0/. should be called 0/1 for heterozygous insertion
                elif(geno_list[0] == "3.0/."):
                    geno_list[0] = "0/1"
                # 3.5/. and 4.0 should both be called 1/1 for homozygous insertion
                elif(geno_list[0] == "3.5/." or geno_list[0] == "4.0/."):
                    geno_list[0] = "1/1"
            line_list[a] = ":".join(geno_list)
                
        # prepare output
        output_string = "\t".join(line_list)
        output_string = output_string.strip("\n")
        output_string = output_string + "\n"
        file_out.write(output_string)
            
    
    file_in.close()
    file_out.close()
    
#vcf_manipulate("../Data/CNV.vcf", "../Data/CNV2.vcf")
tsv_manipulate("../Data/CNV2_anno.tsv","../Data/CNV2_anno_corrected.tsv")
