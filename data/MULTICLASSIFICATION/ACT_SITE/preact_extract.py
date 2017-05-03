#!/usr/bin/env python
import sys
import re
import urllib
from Bio import SeqIO
filename=sys.argv[1]
pep_list=[]
for line in open(filename).readlines():
	line=line.strip()
	name=line.split()[0]
	cnum=line.split()[1]
	peptide=line.split()[2][5:36]
	if peptide[15]=='C':
		if peptide not in pep_list:
			pep_list.append(peptide)
			print name,cnum,peptide
