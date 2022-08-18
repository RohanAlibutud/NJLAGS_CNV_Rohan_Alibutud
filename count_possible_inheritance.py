import csv
import collections
import itertools
input_path = "../Data/sample_summary.tsv"



ped_dict = dict()
count_dict = dict()
fam_summary = collections.defaultdict(dict)

pheno_list = ["ASD","ADHD","LI","RI","SRS"]
pattern_list = ["dom","rec"]
for pattern in pattern_list:
    count_dict[pattern] = {}
    for pheno in pheno_list:
        count_dict[pattern][pheno] = set()

#Read in ped file
with open(input_path, 'r') as read_obj:
    reader = csv.DictReader(read_obj,delimiter='\t')
    #fid    iid fatid   matid   sex ASD LI  RI  SRS ADHD

    for row in reader:
        famid = row["Family ID"]
        iid = row["PARTID"]
        fatid = row["Father ID"]
        matid = row["Mother ID"]
        ASD = row["ASD"]
        LI = row["LI"]
        RI = row["RI"]
        SRS = row["SRS Cat"]
        ADHD = row["ADHD_evidence"]
        if iid not in ped_dict:
            ped_dict[iid] = {"family": famid,"father": fatid, "mother":matid, "ASD":ASD,"LI":LI,"RI":RI,"SRS":SRS,"ADHD":ADHD}
        if famid not in fam_summary:
            for pheno in pheno_list:
                fam_summary[famid][pheno] = 0
            fam_summary[famid]["trio"] = "N"
            fam_summary[famid]["samples"] = 0


for iid in ped_dict:
    famid = ped_dict[iid]["family"]
    fam_summary[famid]["samples"] += 1
    for pheno in pheno_list:
        if ped_dict[iid][pheno] == "2":
            fam_summary[famid][pheno] += 1
    if ped_dict[iid]["father"] != "0" and ped_dict[iid]["mother"] != "0":
        mid = ped_dict[iid]["father"]
        fid = ped_dict[iid]["mother"]
        famid = ped_dict[iid]["family"]


        if mid in ped_dict and fid in ped_dict:
            if fam_summary[famid]["trio"] == "Y":
                fam_summary[famid]["trio"] = "N"
            else:
                fam_summary[famid]["trio"] = "Y"
  
            for pheno in pheno_list:
                if ped_dict[iid][pheno] == "2":                    
                    if ped_dict[mid][pheno] == "2" or ped_dict[fid][pheno] == "2":
                            count_dict["dom"][pheno].add(famid)
                            #fam_assigned_dom.add(famid)
                            name = 'dom '+ pheno
                            fam_summary[famid][name] = 1
                    else:
                        count_dict["rec"][pheno].add(famid)
                        name = 'rec '+ pheno
                        fam_summary[famid][name] = 1
                        #fam_assigned_rec.add(famid)


for pattern in pattern_list:
    for pheno in pheno_list:
        size = len(count_dict[pattern][pheno])
        print("%s-%s %d" % (pattern,pheno,size))

fields = ['famid','ASD','ADHD','LI','RI','SRS','samples','trio','dom ASD', 'dom ADHD', 'dom LI', 'dom RI', 'dom SRS','rec ASD', 'rec ADHD', 'rec LI', 'rec RI', 'rec SRS']
with open('../Data/family_summary.csv','w',newline='') as f:
    w = csv.DictWriter(f,fields)
    w.writeheader()
    for key,value in fam_summary.items():
        row = {'famid': key}
        row.update(value)
        w.writerow(row)

#print(len(fam_assigned_dom))
#print(len(fam_assigned_rec))

    
        #if row['AF_3'] == 'None' and row['AF_2'] == 'None':
        #    if float(row['sampleAF']) <= 0.001:
        #        writer.writerow(row)

        #elif (row['AF_2'] != 'None' and float(row['AF_2']) <= 0.001) or (row['AF_3'] != 'None' and float(row['AF_3']) <= 0.001):
        #    writer.writerow(row)
