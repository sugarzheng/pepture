#!/usr/bin/env python
import os, random
import numpy as np
import h5py
import math
import sys
from Bio import SeqIO
AA_abbre={'D':'aspartate','E':'glutamate','C':'cysteine','H':'histidine','S':'serine'}

def extract_peptide(fulseq,name,num):
	cnum=int(num)
	low_=max(cnum-21,0)
	high_=min(len(fulseq),cnum+20)
	line=name+"\t"+str(cnum)+"\t"+(low_-(cnum-21))*'-'+fulseq[low_:high_]+(cnum+20-high_)*'-'+'\n'
	return line

def download_fasta(database):
	if os.path.exists(database):
		pass
	else:
		if not os.path.exists(database+".gz"):
			os.system("wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz -P data/")
		print("unziping data/uniprot_sprot.fasta.gz")
		os.system("gunzip data/uniprot_sprot.fasta.gz ")	


def extract_active(AA):
	AA_name=AA_abbre[AA]
	print("start extracting "+AA_name+" active site from active_site.dat file")
	pep_list=[]
	filename="data/active_site.dat"
	output=open("data/"+AA_name+"/active.dat",'w')
	act_list=open(filename).readlines()
	database='data/uniprot_sprot.fasta'
	donwload_fasta(database)
	records = SeqIO.parse(database,"fasta")
	record_dict=SeqIO.to_dict(records)
	print("the fasta dictionary has been prepared")
	for line in act_list:
		line=line.strip()
		gene_name="sp|"+line.split(";")[0]+"|"+line.split()[1]
		num=line.split()[2]
		try:
			fulseq=str(record_dict[gene_name].seq)
		except:
			continue
		if fulseq[int(num)-1]==AA:
			write_line=extract_peptide(fulseq,gene_name,num)
			peptide=write_line.strip().split()[2][5:36]
			if peptide not in pep_list:
				pep_list.append(peptide)
				output.write(write_line)
	output.close()
	print(AA_name+" active site has been extracted.")
	print("start extracting not active "+AA_name+" site.")
	not_output=open("data/"+AA_name+"/not_active.dat",'w')
	special_act_list=open("data/"+AA_name+"/active.dat").readlines()
	not_active_list=[]
	for line in special_act_list:
		line=line.strip()
		gene_name=line.split()[0]
		num=line.split()[1]
		try:
			fulseq=str(record_dict[gene_name].seq)
		except:
			print("bad data "+gene_name)
			continue
		for (i,aa) in enumerate(fulseq):
			if aa == AA:
				pep=extract_peptide(fulseq,gene_name,i+1)
				if pep not in special_act_list:
					if pep.strip().split()[2][10:30] in not_active_list:
						pass
					else:
						not_active_list.insert(0,pep.strip().split()[2][10:30])
						not_output.write(pep)
	not_output.close()	
	print(AA_name+"'s not active site has been extracted.")

def cut_line(line,num):
	line=line.strip()
	if len(line.split()[2])-1<num:
		print('the maxium length is '+str(len(line.split()[2]))+" you can not use the length number "+str(num))
		sys.exit()
	peptide_cut=(len(line.split()[2])-1-num)/2
	line=line.split()[0]+"\t"+line.split()[1]+"\t"+line.split()[2][peptide_cut:(len(line.split()[2])-peptide_cut)]
	return line

def data_append(list_file,(train,test,validation),label):
	for (i,line) in enumerate(list_file):
		if '-' not in line.split()[2]:
			if i%3 != 0:
				train.append(line.strip()+" "+str(label))
			else:
				if i%6!=0:
					test.append(line.strip()+" "+str(label))
				else:
					validation.append(line.strip()+" "+str(label))
	return (train,test,validation)

def data_append_num(list_file,(train,test,validation),num,label):
	list_length=len(list_file)
	for (i,line) in enumerate(list_file):
		line=cut_line(line,num)
		if '-' not in line.split()[2]:
#			if line not in (train + test + validation):
			if i%3 != 0:
				train.append(line.strip()+" "+str(label))
			else:
				if i%6!=0:
					test.append(line.strip()+" "+str(label))
				else:
					validation.append(line.strip()+" "+str(label))
	return (train,test,validation)

def reverse_add(dataset):
	reverse_list=list()
	for line in dataset:
		line=line.strip()
		line=line.split()[0]+"_reverse"+" "+line.split()[1]+" "+line.split()[2][::-1]
		reverse_list.append(line)
	return reverse_list
	

###########################################start inputing data
#def loadmulticlassification():
#	active_site = open("data/MULTICLASSIFICATION/ACT_SITE/input.dat", "r").readlines()
#	COMBINE=open("data/MULTICLASSIFICATION/combine/input.dat",'r').readlines()
#	disulfide=open("data/MULTICLASSIFICATION/INNERDISULFIDE/input.dat",'r').readlines()
#	train=[]
#	test=[]
#	validation=[]
#	random.shuffle(active_site)
#	random.shuffle(COMBINE)
#	random.shuffle(disulfide)
#	(train,test,validation)=data_append(COMBINE,(train,test,validation),0)
#	(train,test,validation)=data_append(active_site,(train,test,validation),1)
#	(train,test,validation)=data_append(disulfide,(train,test,validation),2)
#	random.shuffle(train)
#	random.shuffle(test)
#	random.shuffle(validation)
#	return (train,validation,test,3)

def loadtestdata():
	active_site = open("data/small_test_data/ACT_SITE/input.dat", "r").readlines()
	COMBINE=open("data/small_test_data/combine/input.dat",'r').readlines()
	#disulfide=open("data/small_test_data/INNERDISULFIDE/input.dat",'r').readlines()
	train=[]
	test=[]
	validation=[]
	(train,test,validation)=data_append(COMBINE,(train,test,validation),0)
	(train,test,validation)=data_append(reverse_add(active_site),(train,test,validation),0)
	(train,test,validation)=data_append(reverse_add(COMBINE),(train,test,validation),0)
	(train,test,validation)=data_append(active_site,(train,test,validation),1)
	#(train,test,validation)=data_append(disulfide,(train,test,validation),2)
	return (train,validation,test,2)

#def loadbinaryclassification():
#	active_site = open("data/ACTIVE_OR_NOT/active.dat", "r").readlines()
#	not_active_site=open("data/ACTIVE_OR_NOT/not_active.dat", "r").readlines()
#	###########reverse sequence as negative control
#	not_active_site=reverse_add(active_site,not_active_site)
#	print(not_active_site)
#	train=[]
#	test=[]
#	validation=[]
#	random.shuffle(active_site)
#	random.shuffle(not_active_site)
#	(train,test,validation)=data_append(not_active_site,(train,test,validation),0)
#	(train,test,validation)=data_append(active_site,(train,test,validation),1)
#	random.shuffle(train)
#	random.shuffle(test)
#	random.shuffle(validation)
#	return(train,validation,test,2)
#
def length_classification(file1,file2,num):
	###num means the length of peptide except the center residue,num must be an even number
	active_site = open(file1, "r").readlines()
	not_active_site=open(file2, "r").readlines()
	train=[]
	test=[]
	validation=[]
	random.shuffle(active_site)
	random.shuffle(not_active_site)
	(train,test,validation)=data_append_num(not_active_site,(train,test,validation),num,0)
	(train,test,validation)=data_append_num(reverse_add(not_active_site),(train,test,validation),num,0)
	(train,test,validation)=data_append_num(reverse_add(active_site),(train,test,validation),num,0)
	(train,test,validation)=data_append_num(active_site,(train,test,validation),num,1)
	random.shuffle(train)
	random.shuffle(test)
	random.shuffle(validation)
	return(train,validation,test,2)

def download_data(AA,num):
	AA_name=AA_abbre[AA]
	uniprot_file="data/uniprot_sprot.dat"
	if os.path.exists(uniprot_file):
		pass
	else:
		if not os.path.exists(uniprot_file+".gz"):
			os.system("wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz -P data/")
		print("unzipping uniprot_sprot.dat.gz")
		os.system("gunzip data/uniprot_sprot.dat.gz ")
	if os.path.exists("data/active_site.dat"):
		pass
	else:
		print("preparing for all the active_site dataset in the uniprot_sprot database")
		os.system("awk \'{if($1==\"ID\")ID=$2;if($1==\"AC\")AC=$2;if($1==\"FT\" && $2==\"ACT_SITE\")print AC,ID,$3}\' data/uniprot_sprot.dat >data/active_site.dat")
	AA_file_1="data/"+AA_name+"/active.dat"
	AA_file_0="data/"+AA_name+"/not_active.dat"
	if os.path.exists(AA_file_1) and os.path.exists(AA_file_1) and os.path.getsize(AA_file_1) and os.path.getsize(AA_file_0):
		print("the "+AA_name+" file has already prepared")
	else:
		if os.path.exists("data/"+AA_name):
			pass
		else:
			os.mkdir("data/"+AA_name)
		extract_active(AA)
		
def load_data(AA,num):
	download_data(AA,num)
	AA_name=AA_abbre[AA]
	return length_classification('data/'+AA_name+'/active.dat','data/'+AA_name+'/not_active.dat',num)	
