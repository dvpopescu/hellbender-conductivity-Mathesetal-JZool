# hellbender-conductivity-Mathesetal-JZool
code and data for study on effects of conductivity on larval eastern hellbenders

This code implements analyses described in Mathes, S., S. Kuchta, K. Johnson and V. D. Popescu (2026) Developmental, physiological and fitness consequences of chronic exposure to high conductivity in larval Eastern Hellbenders. Journal of Zoology 

There are 2 scripts:
1) "survival analysis.R" - script performing a survival analysis using the file "surv.csv"
2) "mort_morph_perf.R" - script evaluating differences in survival (based on count data), locomotion performance and multiple metrics of morphology between low and high conductivity treatments and between conditions before and after switching animals from low-to-high and high-to-low conductivity treatments. The data used in this script are: "mort_Oct23.csv" (mortality data at the end of experiement and at switching), "morph_Oct23.csv" (morphology data at the end of experiement and at switching), "PerfNov2.csv" (locomotion performance data at the end of experiement and at switching)    

There are 2 additional Excel files that present the data and results of the CORT analysis and the HNE (lipid peroxidation) analysis: "BenderCORT.xls" and "BenderHNE.xls"
