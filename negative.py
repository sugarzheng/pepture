#!/usr/bin/env python
import sys
import urllib
import re
from Bio import SeqIO

def extract_peptide(fulseq,name,num):
	cnum=int(num)
	low_=max(cnum-21,0)
	high_=min(len(fulseq),cnum+20)
	temp=name+"\t"+str(cnum)+"\t"+(low_-(cnum-21))*'-'+fulseq[low_:high_]+(cnum+20-high_)*'-'+'\n'
	print temp
	return temp


filename=sys.argv[1]
AA=sys.argv[2]
output=open("negative_aspartate.dat",'w')
act_list=open(filename).readlines()
database='/home/whbpt/database/uniprot_sprot.fasta'
records = SeqIO.parse(database,"fasta")
record_dict=SeqIO.to_dict(records)
pep_list=[]
for line in act_list:
	line=line.strip()
	gene_name=line.split()[0]
	num=line.split()[1]
	try:
		fulseq=str(record_dict[gene_name].seq)
		count=0		
		for (i,aa) in enumerate(fulseq):
			if aa == AA:
				pep=extract_peptide(fulseq,gene_name,i+1)
				if pep not in act_list and 'X' not in pep and 'U' not in pep:
					if pep[10:30] not in pep_list:
						pep_list.append(pep[10:30])
						output.write(pep)
	except:
		continue
output.close()

