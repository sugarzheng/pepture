#!/usr/bin/env python
from __future__ import print_function
import os
import h5py
from loaddata import *
from peptotensor import *	
def parse_part(AA):
	active_site = open("results/"+AA_abbre[AA]+"/active_site.dat", "r")
#	input_site=open("data/"+AA_abbre[AA]+"/active.dat",'r').readlines()
	input_site=open("data/"+AA_abbre[AA]+"/homology_result_uniq",'r').readlines()
	input_list1=[]
	input_list2=[]
	for line in input_site:
		line=line.strip()
		line=line.split()
		input_list1.append(line[2])
		input_list2.append(line[0].split("|")[2].split("_")[0])
	output=open("results/"+AA_abbre[AA]+"/to_be_tested_site.dat",'w')
	for line in active_site:
		line=line.strip()
		if line.split()[1].split("|")[2].split("_")[0] not in input_list2:
			if line.split()[3] not in input_list1:
				flag=1
				for aa in line.split()[3]:
					if aa not in aa_list:
						flag=0
				if flag==1:
					print(line)
					output.write(line+"\n")
	output.close()
#def homology_part(AA):
#	active_site = "data/"+AA_abbre[AA]+"/active.dat"
#	inputhandle=open(active_site,'r').readlines()
#	homology_path="data/homology_tempfile"
#	if not os.path.exists(homology_path):
#		os.mkdir(homology_path)
#	else:
#		pass
#	for line in inputhandle:
#		line=line.strip()
#		gene=line.split()[0].split("|")[1]
#		num=line.split()[1]
#		peptide=line.split()[2]
#		name=gene+"_"+num
#		filesequence = open(name+".fst", "w")
#		filesequence.write("".join(peptide))
#		filesequence.close()
#		homology_file=homology_path+"/"+ name + "_uniprot.blast"
##		pwd_path=os.popen("pwd").readlines()[0].strip()
#		fasta_path="data/Swiss-Prot/uniprot_sprot"
#		os.chdir(homology_path)
#		os.system("blastp -query " + name + ".fst" + "  -db "+fasta_path+"  -out "+homology_file + " -outfmt 0 -evalue 0.001")
#		uniprotlines=open(homology_file).readlines()	
#		nstart=0
#		for nl, uniprotline in enumerate(uniprotlines):
#			if "> sp" in uniprotline:
#				nstart = 1
#				species = uniprotline.split("_")[1].split()[0]
#				genename=uniprotline.split()[1]
#			if nstart == 1:
#				if "Query" in uniprotline:
#					cnum_start = int(uniprotline.split()[1])
#					cnum_end = int(uniprotline.split()[3])
#				if "Sbjct" in uniprotline:
#					ccnum_start=int(uniprotline.split()[1])
#					ccnum_end=int(uniprotline.split()[3])
#					peptide=uniprotline.split()[2]
#					hyphon_count=peptide[0:21-cnum_start].count('-')	
#					try:
#						cysteine=uniprotline.split()[2][21-cnum_start]
#						cysteine_num=ccnum_start+21-cnum_start- hyphon_count
#						if cysteine==AA:
#							peptide=(cnum_start-1)*"-"+uniprotline.split()[2]+(41-cnum_end)*"-"
#							if 'X' not in peptide and '-' not in peptide:
#								print(genename,cysteine_num,peptide[0:41])#,ccnum_start,ccnum_end
#					except:
#						continue
#		os.chdir('..')
