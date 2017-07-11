# Post-processing for S-PrediXcan and multiple-tissue pipelines 

## Introduction 
+ A pipeline is used to postprocess S-PrediXcan and multiple-tissue pipeline's results 

## Prerequisites
+  [Python 2.7+](http://www.python.org/download/)
+  [R 3.0+](http://www.r-project.org/)
+  [rpy2](http://rpy2.readthedocs.io/en/version_2.7.x/)
+  [annotables](https://github.com/stephenturner/annotables#how)
+  [dplyr](https://github.com/hadley/dplyr)
+  [qqman](https://github.com/stephenturner/qqman)
+  [ggplot2](https://github.com/hadley/ggplot2)

## Installation
```bash 
git clone https://github.com/jiamaozheng/S-PrediXcan_Post-processing_Pipeline
``` 

## Input and LocusZoom Directories    
+ Create a new folder - `input` which includes the following tools/files: 
   * Outputs either from S-PrediXcan or multiple_tissue pipeline (**required**)
   * [Prediction models](http://hakyimlab.org/predictdb/) (**optional**, only for analyzing S-PrediXcan output) 
   * `gwas_snp.txt` file (tab-delimited format) prepared according to [instructions](http://genome.sph.umich.edu/wiki/LocusZoom_Standalone) and [a concrete example](https://s3.amazonaws.com/imlab-jiamaoz/shared/gwas_snp.txt) (**optional**, only for generating locuszoom plots)
   * [plink](http://pngu.mgh.harvard.edu/~purcell/plink/) (**optional**, only for generating locuszoom plots)

+ [Download LocusZoom (37.5GB)](http://genome.sph.umich.edu/wiki/LocusZoom_Standalone)(**optional**, only for generating locuszoom plots)


## Command Line Parameters 
  Argument              |  Abbre  | Required | Default  | Description  
  ----------------------| ------- | -------- | -------- | ------------------------
  --project_name	    |  -p     |   Yes    |  None    | project name (e.g ovarian_cancer)
  --multiple_tissue      |  -m     |   No     |  ''      | Type "true" for multiple_tissue
  --locuszoom      |  -l     |   No     |  ''      | Type "false" for locuszoom

## Running Pipeline  
**Example 1: To analyze outputs from S-PrediXcan pipeline**
 ```bash 
 python run.py -p ovarian_cancer
 ``` 

**Example 2: To analyze outputs from S-PrediXcan pipeline without locuzoom plot analysis**
 ```bash 
 python run.py -p ovarian_cancer -l false  
 ``` 

**Example 3: To analyze outputs from multiple_tissue pipeline**
 ```bash 
 python run.py -p ovarian_cancer -m true 
 ``` 

**Example 4: To analyze outputs from multiple_tissue pipeline without locuzoom plot analysis**
 ```bash 
 python run.py -p ovarian_cancer -m true -l false  
 ``` 

## Outputs 
 + `../log/` - logs 
 + `../out/` - tables and figures 
    * annotation `.csv`
    * top sigificant genes `.csv`
    * top sigificant genes with snps `.csv` 
    * manhattan_plot `.pdf`
    * qq_plot `.pdf`
    * region_plot `.pdf`
    * bubble_plot `.pdf`
    * locuszoom_plot `.pdf`