#!/usr/bin/env python

''' Edit a table generated with snpSift from a snpEff-annotated VCF file of five individuals, including all genotypes printed out in columns 5-9 (1-based) and the depths of each of these genotypes in columns 10-14. 
Genotypes of a depth <3 are set to missing.

Usage: python 1_filter_table_min_depth3_set_missing.py table.txt table_dp3.txt

'''

import sys

fileIn=open(sys.argv[1], 'r')
fileOut=open(sys.argv[2], 'w')

lines = fileIn.readlines()

for l in lines:
	if l.startswith('#'):
		print >> fileOut, l.rstrip('\n')
	else:
		col4=l.split()[:4] #split lines before and after columns containing genotypes and depths and save in lists (to be merged in the end with edited genotype list)
		gt=l.split()[4:9]
		dpStr=l.split()[9:14]
		col14=l.split()[14:]
		dpInt=[int(i) for i in dpStr] #convert list items in depth list from string to integer
		gtIndex=[] #collect indices of genotypes of depth < 3
		for i, j in enumerate(dpInt):
			if j < 3:
				gtIndex.append(i) #add index of depth < 3 to list gtIndex
				for k in gtIndex:
					gt[k] = './.' #access genotype in list 'gt' with index 'k' and set it to './.'
		newCol4 = '\t'.join(col4).strip()
		newGT = '\t'.join(gt).strip()
		newDP = '\t'.join(dpStr).strip()
		newCol14 = '\t'.join(col14).strip()
		print >> fileOut, newCol4+'\t'+newGT+'\t'+newDP+'\t'+newCol14 #merge the original line pieces with the replaced section and write to output. 

fileIn.close()
fileOut.close()