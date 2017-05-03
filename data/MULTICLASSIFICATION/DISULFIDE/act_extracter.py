#!/usr/bin/env python
import sys
import urllib
from Bio import SeqIO

def extract_peptide(fulseq,name,num):
	cnum=int(num)
	low_=max(cnum-21,0)
	high_=min(len(fulseq),cnum+20)
	temp=name.split("|")[2]+"\t"+str(cnum)+"\t"+(low_-(cnum-21))*'-'+fulseq[low_:high_]+(cnum+20-high_)*'-'+'\n'
	print temp
	return temp

filename="disulfide.dat"
output=open("result.dat",'w')
act_list=open(filename).readlines()
database='/home/whbpt/database/uniprot_sprot.fasta'
records = SeqIO.parse(database,"fasta")
record_dict=SeqIO.to_dict(records)
for line in act_list:
	line=line.strip()
	gene_name="sp|"+line.split(";")[0]+"|"+line.split()[1]
	num=line.split()[2]
	try:
		fulseq=str(record_dict[gene_name].seq)
		output.write(extract_peptide(fulseq,gene_name,num))
	except:
		continue
output.close()
#records = list(SeqIO.parse(database,"fasta"))
#count=0
#gene_name="sp|"+act_list[count].strip().split(";")[0]+"|"+act_list[count].strip().split()[1]
#for record in records:
#	tag=str(record.id)	
#	print gene_name,tag
#	fulseq=str(record.seq)
#	while gene_name in tag:
#		print "bingo"
#		output.write(extract_peptide(fulseq,gene_name,act_list[count].strip().split()[2]))
#		count=count + 1
#		gene_name=act_list[count].strip().split(";")[0]
#output.close()
#pool= ThreadPool(4)
#result=pool.map(extract_peptide,act_list)	
#pool.close()
#pool.join()
#result_list=list()
#for peptide in dict_pep.keys():
#	gene=dict_pep[peptide].split()[0]
#	try:
#		fulseq=str(record_dict[gene].seq)
#	except:
#		continue
#	peptide=peptide.replace('.',"").replace("-","")
#	cdist=peptide.find('*',0)
#	peptide=peptide.replace('*',"")
#	cnum=fulseq.find(peptide,0)+cdist
#	low_=max(cnum-21,0)
#	high_=min(len(fulseq),cnum+20)
#	temp=gene+"\t"+str(cnum)+"\t"+(low_-(cnum-21))*'-'+fulseq[low_:high_]+(cnum+20-high_)*'-'
#	if temp not in result_list:
#		result_list.append(temp)
#outputhandle=open(filename+".peptide",'w')
#for line in result_list:
#	outputhandle.write(line+"\n")
#outputhandle.close()
