#!/usr/bin/env python
from __future__ import print_function
import os
import h5py
from loaddata import *
from peptotensor import *	
def parse_part(AA):
	active_site = open("results/"+AA_abbre[AA]+"/active_site.dat", "r")
	input_site=open("data/"+AA_abbre[AA]+"/active.dat",'r').readlines()
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
