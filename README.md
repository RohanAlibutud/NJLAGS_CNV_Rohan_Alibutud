NJLAGS_CNV_Rohan_Alibutud
Repository for code used in the analysis of copy number variants in the NJLAGS cohort as part of a project in the Xing Lab of Genomics at Rutgers University

CNV_RunThrough.sh

NJLAGS CNV Prioritization Pipeline

Files needed:
-cnv_call_files: all the starting files taken from Vaidhy's analysis, pre-merges
-NJLAGS_CNV.ped: pedigree file from affected families, now including LI, RI, and SRS as phenotypes
-human_g1k_v37.fasta: FASTA file for build 37 of the human reference genome
-human_g1k_v37.fasta.fai: history file for the above FASTA file
-individuals.txt: text file containing all the PARTIDs for the individuals in the sample

Scripts needed:
-CNV_Scripts.py
-CNV_Merge4.py
-VCFizer2.py
-Manipulator.py
-Geminesque_2022.py
-GMN_Prioritizer.py
-Build_StrVCTVRE.py
-StrVCTVRE.py
-Add_StrVCTVRE.py
-GMN_Summarizer.py

Positions on Alu server
Full pipeline run: /lab01/Projects/Rohan_Projects/CNV_Project/Total_Pathway/Scripts/CNV_RunThrough.sh
Resulting files: /lab01/Projects/Rohan_Projects/CNV_Project/2022/Total_Pathway/Results
