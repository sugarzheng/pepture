#!/usr/bin/env python
predict = open("results/cysteine/predict_num", "r").readlines()
align = open("results/cysteine/input.dat", "r").readlines()
active_site = open("results/cysteine/active_site.dat", "w")
for i in range(len(predict)):
    if predict[i].strip().split()[0] == "1":
		active_site.write(predict[i].strip().split()[1]+"\t"+align[i])
