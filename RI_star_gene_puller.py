gene_list = ["LOC105377785", "LOC101927815", "LOC101927354", "MTRNR2L9", "OR2J2", "LINC02241", "KCNIP4", "TANC1", "DAPL1", "PSG7", "PSG6", "C2CD4A", "C2CD4B", "GOLGA2P11", "LOC101928907", "LOC107984784", "VPS13C", "LINC01619", "GRIP1", "HELB", "KCNA6", "IPP", "LOC105378614", "MICOS10", "MICOS10-NBL1", "NBL1", "RPS14P3", "CAPZB", "HNRNPCL2", "HNRNPCL3", "HNRNPCL4"]
file_in = open("../Data/dispensable_genes_and_muc.txt", "r")
file_list = file_in.read().split("\n")

for i in range(0, len(gene_list)):
	if gene_list[i] in file_list or "LOC" in gene_list[i] or "LINC" in gene_list[i]:
		continue
	else:
		print(gene_list[i])
		