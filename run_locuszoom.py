#!/usr/bin/env python
import os
import sys

# Find locuszoom binary. 
sys.argv[0] = os.path.abspath(sys.argv[0]);
lzbin = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]),"../bin/locuszoom"));

# Run a quick example from the Kathiresan et al. data. 
# cmd = "%(bin)s --metal Kathiresan_2009_HDL.txt --refsnp rs174546 --flank 200kb --build hg19 --source 1000G_March2012 --pop EUR title='My region' geneFontSize=1.1 recombColor='gray'" % {'bin' : lzbin};
cmd = "%(bin)s --metal gwas_snp.txt --hitspec batch_locuszoom.txt --build hg19 --source 1000G_March2012 --pop EUR " % {'bin' : lzbin};
print ("Running: %s" % cmd);
os.system(cmd); 
