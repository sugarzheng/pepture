#!/usr/bin/env python
from __future__ import print_function
import keras
from keras import backend as K
import os
import h5py
from keras.models import load_model
from peptotensor import *
from loaddata import *
#def peptide_extracter(fulseq,cnum):
#	low_=max(cnum-20,0)
#	high_=min(len(fulseq),cnum+21)
#	return (low_-(cnum-20))*'-'+fulseq[low_:high_]+(cnum+21-high_)*'-'
def extract_input(fasta,AA):
	result_path="results/"+AA_abbre[AA]+"/input.dat"
	if os.path.exists(result_path):
		pass
	else:
		if not os.path.exists("results/"+AA_abbre[AA]):
			os.mkdir("results/"+AA_abbre[AA])
		record = SeqIO.parse(open(fasta),'fasta')
		record_dict=SeqIO.to_dict(record)
		outputhandle=open("results/"+AA_abbre[AA]+"/input.dat",'w')
		for gene in record_dict.keys():
			fulseq=str(record_dict[gene].seq)
			for (i,aa) in enumerate(fulseq):
				if aa == AA:
					outputhandle.write(extract_peptide(fulseq,gene,i+1))
		outputhandle.close()
	return open(result_path, "r").readlines()	
	
def predict_part(AA,fasta,peptide_num):
	download_fasta(fasta)
	predict_input = extract_input(fasta,AA) 
	prediction=[]
	num_path="results/"+AA_abbre[AA]+"/predict_num"
	if os.path.exists(num_path):
		os.remove(num_path)
	batch_size = 2048
	img_rows, img_cols = 20, peptide_num+1
	model = load_model('models/'+AA_abbre[AA]+'_test.h5')
	count=0
	for (i,line) in enumerate(predict_input):
		line=cut_line(line,peptide_num)
		if i%batch_size!=0 and i!=(len(predict_input)-1) :	
			prediction.append(line+" 0")
		else:
			count=count+batch_size
			print(count)
			prediction.append(line+" 0")
			(x_prediction, y_prediction)=peptoblosum(prediction)
			x_prediction = x_prediction.reshape(x_prediction.shape[0], img_rows, img_cols, 1)
			classes =model.predict(x_prediction,batch_size=batch_size)
			output=open(num_path,'a')
			for i,num in enumerate(classes.argmax(axis=1)):
				output.write(str(num)+"\t"+str(classes[i][num])+"\n")
			prediction=[]
	output.close()
	predict = open(num_path, "r").readlines()
	assert(len(predict)==len(predict_input))
	active_site = open("results/"+AA_abbre[AA]+"/active_site.dat", "w")
	for i in range(len(predict)):
		if predict[i].strip().split()[0] == "1":
 			active_site.write(predict[i].strip().split()[1]+"\t"+predict_input[i])
	active_site.close()
