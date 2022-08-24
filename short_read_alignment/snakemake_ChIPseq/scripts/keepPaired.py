#!/usr/bin/env python3

# Python 3 script adapted from Python 2 script by Devon Ryan ( https://www.biostars.org/p/95929/ )
# for retaining filtered alignments that consist of both reads in a pair
# It achieves this by checking for matching consecutive read names corresponding to a read pair
  
import csv
import sys

#f = csv.reader(sys.stdin, dialect="excel-tab")
#of = csv.writer(sys.stdout, dialect="excel-tab")
#with open('WT_DMC1_V5_Rep1_ChIP_MappedOn_Athaliana_ONT_RaGOO_v2.0_sort.txt', 'r') as file, open('WT_DMC1_V5_Rep1_ChIP_MappedOn_Athaliana_ONT_RaGOO_v2.0_sort_keeppaired.txt', 'w') as outfile:
f = csv.reader(sys.stdin, dialect=csv.excel_tab)
# Samtools v1.10 requires that fields in @PG line not be enclosed in quotation marks
of = csv.writer(sys.stdout, dialect=csv.excel_tab, quoting=csv.QUOTE_NONE, escapechar='"', doublequote=False)
last_read = None
for line in f:
     # Take care of the header
     if(line[0][0] == "@"):
         of.writerow(line)
         continue 
     if(last_read == None):
         last_read = line
     else:
         if(last_read[0] == line[0]):
             of.writerow(last_read)
             of.writerow(line)
             last_read = None
         else:
             last_read = line
