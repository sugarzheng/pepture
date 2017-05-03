#!/usr/bin/env python
import os, random
import numpy as np
import h5py
import math
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
	for (i,line) in enumerate(COMBINE):
		if '-' not in line.split()[2]:
			if i%3 != 0:
				train.append(line.strip()+" 0")
			else:
				if i%6!=0:
					test.append(line.strip()+" 0")
				else:
					validation.append(line.strip()+" 0")
	
	for (i,line) in enumerate(active_site):
		if '-' not in line.split()[2]:
			if i%3 != 0:
				train.append(line.strip()+" 1")
			else:
				if i%6!=0:
					test.append(line.strip()+" 1")
				else:
					validation.append(line.strip()+" 1")
	for (i,line) in enumerate(disulfide):
		if '-' not in line.split()[2]:
			if i%3 != 0:
				train.append(line.strip()+" 2")
			else:
				if i%6!=0:
					test.append(line.strip()+" 2")
				else:
					validation.append(line.strip()+" 2")
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
	random.shuffle(active_site)
	random.shuffle(COMBINE)
	random.shuffle(disulfide)
	for (i,line) in enumerate(COMBINE):
		if '-' not in line.split()[2]:
			if i%3 != 0:
				train.append(line.strip()+" 0")
			else:
				if i%6!=0:
					test.append(line.strip()+" 0")
				else:
					validation.append(line.strip()+" 0")
	
	for (i,line) in enumerate(active_site):
		if '-' not in line.split()[2]:
			if i%3 != 0:
				train.append(line.strip()+" 1")
			else:
				if i%6!=0:
					test.append(line.strip()+" 1")
				else:
					validation.append(line.strip()+" 1")
	for (i,line) in enumerate(disulfide):
		if '-' not in line.split()[2]:
			if i%3 != 0:
				train.append(line.strip()+" 2")
			else:
				if i%6!=0:
					test.append(line.strip()+" 2")
				else:
					validation.append(line.strip()+" 2")
	random.shuffle(train)
	random.shuffle(test)
	random.shuffle(validation)
	return (train,validation,test,3)

def loadbinaryclassification():
	active_site = open("data/ACTIVE_OR_NOT/active.dat", "r").readlines()
	non_active_site=open("data/ACTIVE_OR_NOT/not_active.dat", "r").readlines()
	train=[]
	test=[]
	validation=[]
	random.shuffle(active_site)
	random.shuffle(non_active_site)
	for (i,line) in enumerate(active_site):
		if '-' not in line.split()[2]:
			if i%3 != 0:
				train.append(line.strip()+" 1")
			else:
				if i%6!=0:
					test.append(line.strip()+" 1")
				else:
					validation.append(line.strip()+" 1")
	for (i,line) in enumerate(non_active_site):
		if '-' not in line.split()[2]:
			if i%3 != 0:
				train.append(line.strip()+" 0")
			else:
				if i%6!=0:
					test.append(line.strip()+" 0")
				else:
					validation.append(line.strip()+" 0")
	random.shuffle(train)
	random.shuffle(test)
	random.shuffle(validation)
	return(train,validation,test,2)
