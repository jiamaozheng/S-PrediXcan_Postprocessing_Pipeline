# helper functions

'''
################################
#### Software prerequistes 
################################
'''
# python packages: 
# rpy2:           pip install rpy2

# R packages: 
# annotables:     install.packages("devtools")
#                 devtools::install_github("stephenturner/annotables")
# dplyr:          install.packages("dplyr")
# ggplot2:        install.packages("ggplot2") 
# qqman:          install.packages("qqman") 

# Load modules and libraries 
import os, sqlite3, time, sys, pandas, glob, shutil, subprocess, logging
from datetime import datetime
from rpy2.robjects import r
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2 import robjects
import uuid as myuuid 
pandas2ri.activate()
importr('annotables')
qqman = importr('qqman')
importr('ggplot2')
importr('dplyr',  on_conflict="warn")

# get Log
def getLog(logger, log_path):
    logger = logging.getLogger()
    fhandler = logging.FileHandler(filename=log_path, mode='w')
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fhandler.setFormatter(formatter)
    logger.addHandler(fhandler)
    logger.setLevel(logging.INFO)

# Pretty string for a given number of seconds.
def timeString(seconds):
  tuple = time.gmtime(seconds);
  days = tuple[2] - 1;
  hours = tuple[3];
  mins = tuple[4];
  secs = tuple[5];
  if sum([days,hours,mins,secs]) == 0:
    return "<1s";
  else:
    string = str(days) + "d";
    string += ":" + str(hours) + "h";
    string += ":" + str(mins) + "m";
    string += ":" + str(secs) + "s";
  return string;




