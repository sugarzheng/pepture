#!/usr/bin/env python
import os, random
import numpy as np
import h5py
import math
import sys
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


###########################################start inputing data
def loadmulticlassification():
	active_site = open("data/MULTICLASSIFICATION/ACT_SITE/input.dat", "r").readlines()
	COMBINE=open("data/MULTICLASSIFICATION/combine/input.dat",'r').readlines()
	disulfide=open("data/MULTICLASSIFICATION/INNERDISULFIDE/input.dat",'r').readlines()
	train=[]
	test=[]
	validation=[]
	random.shuffle(active_site)
	random.shuffle(COMBINE)
	random.shuffle(disulfide)
	(train,test,validation)=data_append(COMBINE,(train,test,validation),0)
	(train,test,validation)=data_append(active_site,(train,test,validation),1)
	(train,test,validation)=data_append(disulfide,(train,test,validation),2)
	random.shuffle(train)
	random.shuffle(test)
	random.shuffle(validation)
	return (train,validation,test,3)
def loadtestdata():
	active_site = open("data/small_test_data/ACT_SITE/input.dat", "r").readlines()
	COMBINE=open("data/small_test_data/combine/input.dat",'r').readlines()
	disulfide=open("data/small_test_data/INNERDISULFIDE/input.dat",'r').readlines()
	train=[]
	test=[]
	validation=[]
	(train,test,validation)=data_append(COMBINE,(train,test,validation),0)
	(train,test,validation)=data_append(active_site,(train,test,validation),1)
	(train,test,validation)=data_append(disulfide,(train,test,validation),2)
	return (train,validation,test,3)

def loadbinaryclassification():
	active_site = open("data/ACTIVE_OR_NOT/active.dat", "r").readlines()
	non_active_site=open("data/ACTIVE_OR_NOT/not_active.dat", "r").readlines()
	train=[]
	test=[]
	validation=[]
	random.shuffle(active_site)
	random.shuffle(non_active_site)
	(train,test,validation)=data_append(non_active_site,(train,test,validation),0)
	(train,test,validation)=data_append(active_site,(train,test,validation),1)
	random.shuffle(train)
	random.shuffle(test)
	random.shuffle(validation)
	return(train,validation,test,2)

def length_classification(file1,file2,num):
	###num means the length of peptide except the center residue,num must be an even number
	active_site = open(file1, "r").readlines()
	non_active_site=open(file2, "r").readlines()
	train=[]
	test=[]
	validation=[]
	random.shuffle(active_site)
	random.shuffle(non_active_site)
	peptide_all=[]
	for (i,line) in enumerate(active_site):
		line=cut_line(line,num)
		if '-' not in line.split()[2]:
			if line.split()[2] not in peptide_all:
				peptide_all.append(line.split()[2])
				if i%3 != 0:
					train.append(line+" 1")
				else:
					if i%6!=0:
						test.append(line+" 1")
					else:
						validation.append(line+" 1")
	peptide_all=[]
	for (i,line) in enumerate(non_active_site):
		line=cut_line(line,num)
		if '-' not in line.split()[2]:
			if line.split()[2] not in peptide_all:
				peptide_all.append(line.split()[2])
				if i%3 != 0:
					train.append(line+" 0")
				else:
					if i%6!=0:
						test.append(line+" 0")
					else:
						validation.append(line+" 0")
	random.shuffle(train)
	random.shuffle(test)
	random.shuffle(validation)
	return(train,validation,test,2)
def load_length_classification(num):
	return length_classification('data/ACTIVE_OR_NOT_40_length/active.dat','data/ACTIVE_OR_NOT_40_length/not_active.dat',num)
def load_ABPP_classification(num):
	return length_classification('data/ABPP/active.dat','data/ABPP/not_active.dat',num)
def loadpredictdata(num):
	predict_site=open("results/input.dat", "r").readlines()
	train=[]
	test=[]
	validation=[]
	for (i,line) in enumerate(active_site):
		line=cut_line(line,num)
		train.append(line+" 0")
		test.append(line+" 0")
	return(train,test)
