# S-PrediXcan Post-processing

## Introduction 
+ A pipeline is used to postprocess S-PrediXcan results 

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
+ Create a new fold - `input` which includes the following tools/files: 
   * [plink](http://pngu.mgh.harvard.edu/~purcell/plink/) (**required**)
   * Outputs either from S-PrediXcan or multiple_tissue pipeline (**required**)
   * [Prediction models](http://hakyimlab.org/predictdb/) (optional, only for analyzing S-PrediXcan output) 
   * `gwas_snp.txt` file (tab-delimited format) prepared according to [instructions](http://genome.sph.umich.edu/wiki/LocusZoom_Standalone) and [a concrete example](https://s3.amazonaws.com/imlab-jiamaoz/shared/gwas_snp.txt) (optional, only for generating locuszoom plots)
+ [Download LocusZoom (37.5GB)](http://genome.sph.umich.edu/wiki/LocusZoom_Standalone)(optional software)


## Command Line Parameters 
  Argument              |  Abbre  | Required | Default  | Description  
  ----------------------| ------- | -------- | -------- | ------------------------
  --project_name	    |  -p     |   Yes    |  None    | project name (e.g ovarian_cancer, breast_cancer)
  --multiple_tissue      |  -m     |   NO     |  ''      | multiple_tissue (Type "true", if you would like to analyze outputs from multiple_tissue pipeline)
  --locuszoom      |  -l     |   NO     |  ''      | locuszoom (Type "false", if you don't want to run locuszoom plot analysis)

## Running Pipeline  
+ **Example 1: to analyze outputs from S-PrediXcan pipeline**
 ```bash 
 python run.py -p <your project title>
 ``` 

+ **Example 2: to analyze outputs from S-PrediXcan pipeline without locuzoom plot analysis**
 ```bash 
 python run.py -p <your project title> -l false  
 ``` 

+ **Example 3: to analyze outputs from multiple_tissue pipeline**
 ```bash 
 python run.py -p <your project title> -m true 
 ``` 

+ **Example 4: to analyze outputs from multiple_tissue pipeline without locuzoom plot analysis**
 ```bash 
 python run.py -p <your project title> -m true -l false  
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