"""
VCFizer.py
> turns a tab-delimited file into a proper .vcf
"""

import pysam

# MAIN METHOD
def VCFize():
    print("VCFize method running...")
    
    # open input and output files
    file_in = open("../Data/all_calls.bed", "r")
    file_out = open("../Data/all_calls_unsorted.vcf", "w")
    
    # declare vcf specification
    file_out.write("##fileformat=VCFv4.3\n")
    file_out.write("##fileDate=20220112\n")
    file_out.write("##INFO=<ID=CNV_TYPE,Number=1,Type=String,Description='Type of CNV'>\n")
    file_out.write("##INFO=<ID=END,Number=1,Type=String,Description='Ending position'>\n")
    
    ALT_list = ["##ALT=<ID=CNV,Description='Copy Number Polymorphism'>",
    "##ALT=<ID=CN9,Description='Copy number allele: 9 copies'>",
    "##ALT=<ID=CN8,Description='Copy number allele: 8 copies'>",
    "##ALT=<ID=CN7,Description='Copy number allele: 7 copies'>",
    "##ALT=<ID=CN6,Description='Copy number allele: 6 copies'>",
    "##ALT=<ID=CN5,Description='Copy number allele: 5 copies'>",
    "##ALT=<ID=CN4,Description='Copy number allele: 4 copies'>",
    "##ALT=<ID=CN3,Description='Copy number allele: 3 copies'>",
    "##ALT=<ID=CN2,Description='Copy number allele: 2 copies'>",
    "##ALT=<ID=CN1,Description='Copy number allele: 1 copies'>",
    "##ALT=<ID=CN0,Description='Copy number allele: 0 copies'>"]
    
    contigs_list = ["##contig=<ID=1,assembly=b37,length=249250621>",
    "##contig=<ID=2,assembly=b37,length=243199373>",
    "##contig=<ID=3,assembly=b37,length=198022430>",
    "##contig=<ID=4,assembly=b37,length=191154276>",
    "##contig=<ID=5,assembly=b37,length=180915260>",
    "##contig=<ID=6,assembly=b37,length=171115067>",
    "##contig=<ID=7,assembly=b37,length=159138663>",
    "##contig=<ID=8,assembly=b37,length=146364022>",
    "##contig=<ID=9,assembly=b37,length=141213431>",
    "##contig=<ID=10,assembly=b37,length=135534747>",
    "##contig=<ID=11,assembly=b37,length=135006516>",
    "##contig=<ID=12,assembly=b37,length=133851895>",
    "##contig=<ID=13,assembly=b37,length=115169878>",
    "##contig=<ID=14,assembly=b37,length=107349540>",
    "##contig=<ID=15,assembly=b37,length=102531392>",
    "##contig=<ID=16,assembly=b37,length=90354753>",
    "##contig=<ID=17,assembly=b37,length=81195210>",
    "##contig=<ID=18,assembly=b37,length=78077248>",
    "##contig=<ID=19,assembly=b37,length=59128983>",
    "##contig=<ID=20,assembly=b37,length=63025520>",
    "##contig=<ID=21,assembly=b37,length=48129895>",
    "##contig=<ID=22,assembly=b37,length=51304566>",
    "##contig=<ID=MT,assembly=b37,length=16569>",
    "##contig=<ID=X,assembly=b37,length=155270560>",
    "##contig=<ID=Y,assembly=b37,length=59373566>",
    "##contig=<ID=hs37d5,assembly=b37,length=35477943>"]
    
    for a in range(0, len(contigs_list)):
        file_out.write(contigs_list[a] + "\n")
    
    file_out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype as total number of copy variants'>\n")
    file_out.write("##FORMAT=<ID=PI,Number=1,Type=Integer,Description='PARTID'>\n")
    file_out.write("##FORMAT=<ID=CA,Number=1,Type=String, Description='Caller'>\n")
    file_out.write("##FORMAT=<ID=CH,Number=1,Type=Integer,Description='Chromosome'>\n")
    file_out.write("##FORMAT=<ID=SA,Number=1,Type=Integer,Description='Start position'>\n")
    file_out.write("##FORMAT=<ID=SO,Number=1,Type=Integer,Description='Stop position'>\n")
    file_out.write("##FORMAT=<ID=QL,Number=1,Type=String,Description='Quality reported by CNV filtering'>\n")
    file_out.write("##FORMAT=<ID=PR,Number=1,Type=Integer,Description='Number of probes'>\n")
    
    for a in range(0, len(ALT_list)):
        file_out.write(ALT_list[a] + "\n")
        
    # separate file into lines
    file_list = file_in.read().split("\n")
    
    # change the header
    header_list = file_list[0].split("\t")
    header_list[0] = "#CHROM"
    header_list[1] = "POS"
    header_list[2] = "ID"
    header_list[3] = "REF"
    header_list = pos_insert("ALT",header_list, 5)
    header_list = pos_insert("QUAL",header_list, 6)
    header_list = pos_insert("FILTER",header_list, 7)
    header_list = pos_insert("INFO",header_list, 8)
    header_list = pos_insert("FORMAT",header_list, 9)
    
    # remove trailing tab
    header_list.pop()
    
    header_string = "\t".join(header_list)
    file_out.write(header_string + "\n")
    
    # iterate through input line by line
    for a in range(1, len(file_list)):

        # separate each line into tab-delimited columns
        line_list = file_list[a].split("\t")
        
        # skip empty lines
        if len(line_list) < 2:
            continue
        
        # find CNV genotype to determine del or dup
        for b in range(0, len(line_list)):
            if(":" in line_list[b]):
                CNV_number = line_list[b][0]
                
        # find information for REF allele
        start_pos = int(line_list[1])
        stop_pos = int(line_list[2])
        ref_allele = get_fasta(start_pos, stop_pos)
        
        # find information for ALT allele
        if(line_list[3] == "DELETION ONLY" or line_list[3] == "SPLIT, DELETIONS"):
            alt_allele = "CNV"
        elif(line_list[3] == "DUPLICATION ONLY" or line_list[3] == "SPLIT, DUPLICATIONS"):
            alt_allele = "CNV"
        else:
            if(CNV_number == "0" or CNV_number == "1"):
                alt_allele = "CNV"
            else:
                alt_allele = "CNV"
        
        # fill in empty columns
        line_list[2] = "chr" + line_list[0] + "_" + str(start_pos) + "_" + str(stop_pos)
        line_list[3] = ref_allele # for REF field
        line_list = pos_insert("<CN" + CNV_number + ">", line_list, 5) # for ALT field
        line_list = pos_insert("NULL", line_list, 6) # for QUAL field
        line_list = pos_insert("PASS", line_list, 7) # for FILTER field
        line_list = pos_insert("SVTYPE=" + alt_allele + ";END=" + str(stop_pos), line_list, 8) # for INFO field
        line_list = pos_insert("GT:PI:CA:CH:SA:SO:QL:PR", line_list, 9) # for FORMAT field
        
        # replace alternative genotypes with VCF-correct format
        for c in range(9, len(line_list)):
            if(":" in line_list[c]):
                geno_list = line_list[c].split(":")
                # translate CNV count into genotype
                if(geno_list[0] == "0/."): # homozygous deletion
                    geno_list[0] = "1/1"
                if(geno_list[0] == "1/."): # heterozygous deletion
                    geno_list[0] = "0/1"
                if(geno_list[0] == "3/."): # heterozygous insertion
                    geno_list[0] = "0/1"
                if(geno_list[0] == "4/."): # homozygous insertion
                    geno_list[0] = "1/1"
                line_list[c] = ":".join(geno_list)
                if("|" in line_list[c]):
                    line_list[c] = line_list[c].split("|")[0]
            elif(line_list[c] == "."):
                line_list[c] = "0/0:0:0:0:0:0:null:0"
        
        
        
        # write out line
        line_string = "\t".join(line_list)
        file_out.write(line_string + "\n")
    
    # close input and output files
    file_in.close()
    file_out.close()
    
    print("VCFizing complete.")
    
# INSERT POSITION SUBROUTINE
def pos_insert(x, n_list, pos):
    return n_list[:pos-1]+[x]+n_list[pos-1:]

# FASTA POSITIONS METHOD
def get_fasta(start_pos, stop_pos):
    genome = pysam.Fastafile("../Data/human_g1k_v37.fasta")
    sequence = genome.fetch("1", start_pos, start_pos+1)
    return(sequence)

# SORT .VCF
def VCF_sort():
    # open files
    file_in = open("../Data/all_calls_unsorted.vcf", "r")
    file_out = open("../Data/CNV.vcf", "w")
    
    file_list = file_in.read().split("\n")
    
    skip_line = False
    for a in range(0, len(file_list)-1):
        # skip irrelevant lines
        if(file_list[a][0] == "#"):
            file_out.write(file_list[a] + "\n")
            continue
        else:    
            if(skip_line == False):
                # determine if a swap is necessary
                curr_list = file_list[a].split("\t")
                next_list = file_list[a+1].split("\t")
                # declare boolean
                must_swap = False
                
                # check if chromosomes match
                if(curr_list[0] == next_list[0]):
                    # check if mismatch is present
                    if(int(curr_list[1]) > int(next_list[1])):
                        must_swap = True
                else:
                    must_swap = False
            
                if(must_swap == True):
                    curr_line = "\t".join(curr_list) + "\n"
                    next_line = "\t".join(next_list) + "\n"
                    file_out.write(next_line)
                    file_out.write(curr_line)
                    skip_line = True
                else:
                    file_out.write(file_list[a] + "\n")
                    skip_line = False
            else:
                skip_line = False
                continue
            
    
VCFize()
VCF_sort()
