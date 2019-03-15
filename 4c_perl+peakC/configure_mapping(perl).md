# 4C mapping (perl) pipeline configure



## Introduction

+ A mapping pipeline that maps and filters 4C data.
+ In a 4C experiment DNA fragments are ligated to you fragment of interest, which are amplified using an inverse PCR. These fragments need to mapped back to reference genome. Repetitive fragments need to remove from analysis. Note that non-covered fragment are also of interest to the analysis, since they signal no interaction at this genomic location.



## Require software

1. bwa (bwasw)
2. bedtools (fastaFromBed)
3. Inline::C



## Fastq QC

+ Whenever you want to process fastq file, the first one to do is check the quality of fastq file.

​       ` fastqc file.fastq.gz -o output_path`

+ If quality is good then you can run the 4C mapping, otherwise you need processed raw fastq data to clean data. (For example : adapter error, bar code error) **(Note : cutadapter or trim_galore can be used to covert raw data to clean data)**.



## 4C mapping analysis pipeline

The analysis pipeline consists of three steps.

#### **1. Creating a fragment map**

+ To create a fragment map for you enzyme combination of choice please **run the generate_fragment_map.pl script.** 

+ For a fragment map for DpnII and CviQI of the human genome you would use the following command:

​      `perl generate_fragment_map.pl ucsc.hg19.fa GATC GTAC fragment_map/ `

+ GATC and GTAC is restriction enzyme cutting site of DpnII and CviQI, respectively.
+ The fragment map will be stored in the directory fragment_map/.



#### **2. Identifying the repetitive  fragments**





 

