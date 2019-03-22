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

  ```shell
  fastqc file.fastq.gz -o output_path
  ```

+ If quality is good then you can run the 4C mapping, otherwise you need processed raw fastq data to clean data. (For example : adapter error, bar code error) **(Note : cutadapter or trim_galore can be used to covert raw data to clean data)**.



## 4C mapping analysis pipeline

The analysis pipeline consists of three steps.

#### **1. Creating a fragment map**

+ To create a fragment map for you enzyme combination of choice please **run the generate_fragment_map.pl script.** 

+ For a fragment map for DpnII and CviQI of the human genome you would use the following command:

  ```shell
  perl generate_fragment_map.pl ucsc.hg19.fa GATC GTAC fragment_map/
  ```

+ GATC and GTAC is restriction enzyme cutting site of DpnII and CviQI, respectively.
+ The fragment map will be stored in the directory fragment_map/. (fragment_map/ is an existing directory)



#### **2. Identifying the repetitive  fragments**

+ This step is to filter the fragment map for repetitive fragment, therefore we should map all the fragment back to the genome then selected them from to test whether they are unique or not.

+ To identify the repetitive fragments please **run the getRepeats.pl script.**

  ```shell
  perl getRepeats.pl fragment_map/ restriction_site 49 ucsc.hg19.fa repeat/ [threads]
  ```

+ The results will be placed in the directory repeat/. （repeat/ is an existing directory. In this directory a new directory will created based on the sequence length option that is provided.)

+ **Note : ** 

  1. **restriction_site : ** restriction site of the 1st restriction enzyme.

  2. **49**  , this is the length of ligated fragment including the restriction site. For every different sequencing length for you 4C experiment, you will need to create new repeat map.

     计算方法：如果测序片段长度是150， 引物长度是20, 实验采用的是 4碱基酶，那么这个连接片段的长度 是 150-20-4=134.

  3. reference : reference genome (ucsc.hg19.fa), should be the same as on fragment map is base.

  4. threads : optional parameter, for the amount of threads used in bwa.



#### 3. Splitting FASTQ and mapping to the genome

+ The preprocessing of the data is now finished and you can start map to you data to the genome. the only thing you need is an index file, which contains the minimal information of your 4C experiment. The structure of this file is as follows:

  | Experiment name | primer sequence    | path to reference genome | restriction enzyme1 | restriction enzyme2 | viewpoint chromosome |
  | --------------- | ------------------ | ------------------------ | ------------------- | ------------------- | -------------------- |
  | IRF5            | TCCTCCAGGCCTTGGATC | path/hg19.fa             | GATC                | GTAC                | chr7                 |

+ **Note :**

  1. The title is not necessary that is in order to understand easily. **The index file only need the second line.**

  2. The reference should also have **bwa index.**
  3. Given the current setup it is not possible to mix restriction enzyme combination or reference genomes. If you have multiple genomes or multiple restriction enzyme combinations please create a separate index file for each one.

+ To process and map your data please run **mapping_pipeline.pl**

  ```shell
  perl mapping_pipeline.pl index_file out_dir fastq threads [fragment_map repeat_file]
  ```

  

​     















 

