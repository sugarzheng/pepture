#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from Bio.SubsMat.MatrixInfo import blosum62 as blosum
import math
############################dealing with peptide to images
amino_acids = {'C':0, 'I':1,  'F':2,  'L':3,  'V':4,  'W':5,  'M':6,  'Y':7, 'A':8,  'G':9,  'H':10,  'T':11,  'S':12, \
'P':13,  'Q':14,  'N':15,  'R':16,  'D':17,  'E':18,  'K':19,  '-':20,'B':20,'X':20,'Z':20,'U':20,'O':20}
aa_list=['C','I','F','L','V','W','M','Y','A','G','H','T','S','P','Q','N','R','D','E','K']
def peptovec(vector):
	result_vector=[]
	ytemp=[]
	for line in vector:
		line=line.strip()
		xtemp=[]
		for aa in line.split()[2]:
			aa_list=[0]*20
			aa_list[int(amino_acids[aa])]=1
			xtemp.append(aa_list)
		ytemp.append(line.split()[-1])
		result_vector.append(np.array(xtemp).T)
	return (np.array(result_vector),np.array(ytemp))
def peptoblosum(vector):
    result_vector=[]
    ytemp=[]
    for line in vector:
        line=line.strip()
        xtemp=[]
        for aa in line.split()[2]:
            aa_list_temp=[0]*20
            for (i,num) in enumerate(aa_list):
                if aa in aa_list:
                    try:
                        aa_list_temp[i]=int((math.exp(float(blosum[(num,aa)]*0.347)))/(math.exp(11*0.347))*255)
                    except:
                        aa_list_temp[i]=int((math.exp(float(blosum[(aa,num)]*0.347)))/(math.exp(11*0.347))*255)
                else:
                    aa_list_temp[i]=int((math.exp(-4*0.347))/(math.exp(11*0.347))*255)
            xtemp.append(aa_list_temp)
        ytemp.append(line.split()[-1])
        result_vector.append(np.array(xtemp).T)
    return (np.array(result_vector),np.array(ytemp))
