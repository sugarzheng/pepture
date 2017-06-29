#!/usr/bin/env python
predict = open("results/LEGPH/predict_num", "r").readlines()
align = open("results/LEGPH/legph_tr.dat", "r").readlines()
#train_data=open("data/LEGPH/active.dat",'r').readlines()
#train_data2=open("data/LEGPH/not_active.dat",'r').readlines()
#train_peptide=[]
#for line in train_data:
#	train_peptide.append(line.strip().split()[2])
#for line in train_data2:
#	train_peptide.append(line.strip().split()[2])
active_site = open("results/LEGPH/active_site.pdat", "w")
for i in range(len(predict)):
    if predict[i].strip().split()[0] == "1":
#		if align[i].strip().split()[2] not in train_peptide:
		active_site.write(predict[i].strip().split()[1]+"\t"+align[i])
