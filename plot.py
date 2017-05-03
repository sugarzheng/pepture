#!/usr/bin/env python
from __future__ import print_function
#import keras
#from keras.datasets import mnist
#from keras.models import Sequential
#from keras.layers import Dense, Dropout, Flatten
#from keras.layers import Conv2D, MaxPooling2D
#from keras import backend as K
import os, random
import numpy as np
from peptotensor import *
from loaddata import *
#from pyheatmap.heatmap import HeatMap
from matplotlib import pyplot as PLT
from matplotlib import cm as CM
from matplotlib import axes
############################dealing with peptide to images
def matrix_map(x_train,name):
	sum_=x_train[0]
#	for line in x_train[1:]:
#		sum_+=line
	train_matrix=open(name,'w')
	for line in sum_[0:20]:
		for ll in line:
			train_matrix.write(str(float(ll))+" ")
			#train_matrix.write(str(float(ll)/float(np.sum(line)))+" ")
		train_matrix.write('\n')
	train_matrix.close()
#(train,validation,test,num_classes)=loadmulticlassification()
(train,validation,test,num_classes)=loadtestdata()
(x_train, y_train)=peptoblosum(train)
(x_test, y_test)=peptoblosum(test)
list_name=['tests/a.dat','tests/b.dat','tests/c.dat']
a_list=[]
b_list=[]
c_list=[]
list_list=[a_list,b_list,c_list]
for i in range(0,len(x_train)):
	for j in range(0,num_classes):
		if int(y_train[i])==j:
			list_list[j].append(x_train[i])
for i in range(0,num_classes):
	matrix_map(list_list[i],list_name[i])
