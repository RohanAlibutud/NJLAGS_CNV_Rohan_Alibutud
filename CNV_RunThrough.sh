#!/bin/bash

#starting files needed:
#axiom_run3_combined_20200306_lbf_v2_bpfix_info.csv
#axiom_run6_20190705_dupfix_all_20200114_bpfix_20200514.csv
#axiom_run7_20200302_gc_all_v2_bpfix_info.csv
#illm_jan2015_run2_combined_20200306_lbf_v2_bpfix_info.csv
#illm_may2017_run3_combined_20200306_lbf_v2_bpfix_info.csv
#psych_may2017_run2_20190703_dupfix_all_20200114_bpfix_20200514.csv
#psych_may2017_run3_20200302_gc_all_v2_bpfix_info.csv
#psych_run3_20190709_dupfix_all_20200114_bpfix_20200514.csv
#psych_run4_20200302_gc_fixed_all_v2_bpfix_info.csv
#qsnp_axiom_run2_20190708_dupfix_all_20200114_bpfix_20200514.csv
#qsnp_psych_jan2015_20190711_dupfix_all_20200114_bpfix_20200514.csv
#qsnp_psych_may2017_run2_20190715_dupfix_all_20200114_bpfix_20200514.csv

#scripts to call:
#CNV_Scripts.py
#CNV_Merge4.py
#VCFizer2.py
#Manipulator.py
#Geminesque_2022.py
#GMN_Prioritizer.py
#GMN_Summarizer.py

#ASSEMBLE STARTING DATA INTO .VCF
echo "Beginning CNV Prioritization pipeline"

#1. merge by caller
echo "Running CNV_Scripts.py"
python CNV_Scripts.py
# output: three sorted batch files: axiom_sorted.tsv, jan_sorted.tsv, may_sorted.tsv
echo "CNV_Scripts.py has completed running"

#2. merge by batch
echo "Running CNV_Merge4.py"
python CNV_Merge4.py
# output: all_merged_4.tsv
echo "CNV_Merge4.py has completed running"

#3. rename file to .bed
mv ../Data/merge_outputs/all_merged_4.tsv ../Data/all_calls.bed
# output: all_calls.bed

#4. reformat to .vcf
echo "Running VCFizer2.py"
python VCFizer2.py
# output: CNV.vcf
echo "VCFizer2.py has completed running"

#5. correct some of the format, make a few changes
echo "Running Manipulator1.py"
python Manipulator1.py
# output: CNV2.vcf
echo "Manipulator1.py has completed running"

#6. annotate .vcf using AnnotSV
export ANNOTSV=/lab01/Tools/AnnotSV
$ANNOTSV/bin/AnnotSV -SVinputFile "../Data/CNV2.vcf" -outputFile "../Data/CNV2_anno"
# output: CNV2_anno.vcf

#7. correct more of the format
echo "Running Manipulator2.py"
python Manipulator2.py
# output = CNV2_anno_corrected.vcf
echo "Manipulator2.py has complete running"

#8. annotate using Geminesque
echo "Running Geminesque_2022.py"
python Geminesque_2022.py
# output: GMN_out_2022_01_14.tsv
echo "Geminesque_2022.py has completed running"

#9. prioritize CNVs
echo "Running GMN_Prioritizer.py"
python GMN_Prioritizer.py
# output: GMN_variants_2022_01_14.tsv
echo "GMN_Prioritizer.py has completed running"

#10. add StrVCTVRE annotation
echo "Running StrVCTVRE shell script"
conda env create -f environment_py2.7.yml
conda activate StrVCTVRE_py_2.7
bash StrVCVTRE_runthrough.sh
echo "StrVCTVRE_runthrough.sh has completed running"

