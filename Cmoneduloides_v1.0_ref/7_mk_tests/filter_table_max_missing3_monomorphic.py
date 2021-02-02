#!/usr/bin/env python

''' Edit a table generated with snpSift from a snpEff-annotated VCF file of five individuals, including all genotypes printed out in columns 5-9 (1-based) and the depths of each of these genotypes in columns 10-14. 
Lines with more than 3 missing genotypes and lines with monomorphic sites are removed. 

Usage: python 2_filter_table_max_missing3_monomorphic.py table_dp3.txt table_dp3_varSites_maxMiss3.txt

'''

import sys

fileIn=open(sys.argv[1], 'r')
fileOut=open(sys.argv[2], 'w')

lines = fileIn.readlines()

for l in lines:
	if l.startswith('#'):
		print >> fileOut, l.rstrip('\n')
	else:
		gt=l.replace('/', ' ') #replace '/' from genotypes by whitespace to be able to save every allele into the list 'al'
		al=gt.split()[4:14] #split original columns 5-9 (genotypes) and save elements in a list
		if '.' in al:
			m = l.count('./.') #count number of missing genotypes per line
			al[:] = [j for j in al if j != '.'] #consider only entries that are not missing for discovery of monomorphic sites
			if len(set(al)) != 1 and m <= 3: #print line only if there are more than two different items 'j' in the list 'al', if all genotypes have a depth of at least 3 and if max. 3 genotypes are missing
				print >> fileOut, l.rstrip('\n')
		else:
			if len(set(al)) != 1: #print line only if there are more than two different items 'j' in the list 'al' and if all genotypes have a depth of at least 3
				print >> fileOut, l.rstrip('\n')

fileIn.close()
fileOut.close()
