"""
Build_StrVCTVRE.py
> takes an annotated .tsv file and extracts positions
> composes positions into .bed file
> inputs .bed file into StrVCTVRE 
"""

# BUILD .BED FILE METHOD
def build_bed(var_input, bed_output):
    
    # open files
    var_in = open(var_input, "r")
    bed_out = open(bed_output, "w")
    
    var_list = var_in.read().split("\n")
    for a in range(1, len(var_list)-1):
        line_list = var_list[a].split("\t")
        
        chrom = line_list[0]
        start = line_list[1]
        stop = line_list[2]
        
        # convert SVTYPE to DEL or DUP
        svtype = line_list[3]
        if(svtype == "<CN3>" or svtype == "<CN4>"):
            svtype = "DUP"
        elif(svtype == "<CN0>" or svtype == "<CN1>"):
            svtype = "DEL"
        
        output_list = [chrom, start, stop, svtype]
        output_string = "\t".join(output_list)
        
        bed_out.write(output_string + "\n")
    
    # close files
    var_in.close()
    bed_out.close()
    
build_bed("../Data/prioritized_variants.tsv", "../Data/GMN.bed")