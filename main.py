#!/usr/bin/env python
from __future__ import print_function
import keras
from keras.utils import plot_model
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv2D, MaxPooling2D
from keras import backend as K
import os, random
import numpy as np
import h5py
from sklearn.datasets import make_classification
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, classification_report, confusion_matrix,average_precision_score
from Bio.SubsMat.MatrixInfo import blosum62 as blosum
import math
from peptotensor import * 
from loaddata import *
###########################################start inputing data
#(train,validation,test,num_classes)=loadbinaryclassification()
num=16
img_rows, img_cols = 20, num+1
(train,validation,test,num_classes)=load_length_classification(num)
#(train,validation,test,num_classes)=load_ABPP_classification(num)
#(train,validation,test,num_classes)=load_serine_classification(num)
#(train,validation,test,num_classes)=loadmulticlassification()
#(train,validation,test,num_classes)=loadtestdata()
############start dnn
batch_size = 2048
epochs = 100
(x_train, y_train)=peptoblosum(train)
(x_validation, y_validation)=peptoblosum(validation)
(x_test, y_test)=peptoblosum(test)
#(x_train, y_train)=peptovec(train)
#(x_validation, y_validation)=peptovec(validation)
#(x_test, y_test)=peptovec(test)

print(x_train.shape[0],x_train.shape[1],x_train.shape[2])
x_train = x_train.reshape(x_train.shape[0], img_rows, img_cols, 1)
x_test = x_test.reshape(x_test.shape[0], img_rows, img_cols, 1)
x_validation = x_validation.reshape(x_validation.shape[0], img_rows, img_cols, 1)
input_shape = (img_rows, img_cols, 1)
print('x_train shape:', x_train.shape)
print(x_train.shape[0], 'train samples')
print(x_test.shape[0], 'test samples')
print(x_validation.shape[0], 'validation samples')
print(x_train.shape[0],x_train.shape[1],x_train.shape[2])
# convert class vectors to binary class matrices
y_train = keras.utils.to_categorical(y_train, num_classes)
y_validation = keras.utils.to_categorical(y_validation, num_classes)
y_test = keras.utils.to_categorical(y_test, num_classes)
model = Sequential()
model.add(Conv2D(32, kernel_size=(3, 3),
				activation='relu',
				padding='valid',
				input_shape=input_shape))
#model.add(Conv2D(64, (3, 3), activation='relu'))
#model.add(MaxPooling2D(pool_size=(2, 2)))
#model.add(Dropout(0.25))
model.add(Flatten())
model.add(Dense(512, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(512, activation='relu'))
model.add(Dropout(0.5))
#model.add(Dense(256, activation='relu'))
#model.add(Dropout(0.5))
#model.add(Dense(256, activation='relu'))
#model.add(Dropout(0.5))
#model.add(Dense(256, activation='relu'))
#model.add(Dropout(0.5))
model.add(Dense(num_classes, activation='softmax'))

model.compile(loss=keras.losses.categorical_crossentropy,
			optimizer=keras.optimizers.Adadelta(),
#			optimizer='sgd',
			metrics=['categorical_accuracy'])

model.fit(x_train, y_train, batch_size=batch_size, epochs=epochs,#class_weight=class_weight,
			verbose=1, validation_data=(x_validation, y_validation))
keras.callbacks.TensorBoard(log_dir='./logs')
score = model.evaluate(x_test, y_test, verbose=0)
plot_model(model, to_file='model.png')
print('validation loss:', score[0])
print('validation accuracy:', score[1])
print('validation ?:', score)
model.save('models/cysteine_active_reverse20170626.h5')
accuracy=[]
F1_score=[]
Recall=[]
Precision=[]
weighted_prediction=model.predict(x_validation)
#print(classification_report(weighted_prediction,y_validation))
for i in range(0,num_classes):
	l=y_validation[:,i]
	p=np.array([round(x[i]) for x in weighted_prediction])
	accuracy.append(accuracy_score(l,p))
	F1_score.append(f1_score(l,p))
	Recall.append(recall_score(l,p))
	Precision.append(precision_score(l,p))
print('Accuracy','F1 score','Recall','Precision')
for i in range(0,num_classes):
	print(accuracy[i],F1_score[i],Recall[i],Precision[i])

