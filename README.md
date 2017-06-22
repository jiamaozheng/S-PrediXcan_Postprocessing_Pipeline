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
+ [Download LocusZoom (37.5GB)](http://genome.sph.umich.edu/wiki/LocusZoom_Standalone) 
+ Create a new fold - `input` include the following tools/files:  
   * [Prediction models](http://hakyimlab.org/predictdb/)
   * [plink](http://pngu.mgh.harvard.edu/~purcell/plink/)
   * S-PrediXcan outputs  
   * Prepare `gwas_snp.txt` file (tab-delimited format) as described [here](http://genome.sph.umich.edu/wiki/LocusZoom_Standalone). [A concrete example](https://s3.amazonaws.com/imlab-jiamaoz/shared/gwas_snp.txt) 

## Running Pipeline  
+ Open the terminal, and execute:
 ```python MetaXcanPostprocessing.py <your project title>``` 