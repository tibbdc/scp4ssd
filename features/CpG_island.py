#!/home/zhangjq/anaconda3/envs/september/bin/python
import math
import sys,time
import operator
from copy import deepcopy as dcopy
from functools import reduce
# import pandas as pd
# import numpy as np
# import scipy as sp
log = math.log

OWNZERO = 10**(-30)

training_file = './CpG_training.txt'
result_training = './CpG_training_result.txt'

def read_in(file1):
	observed=' '
	with open(file1) as f:
		observed = observed + reduce(lambda x,y : x.rstrip()+y.rstrip(), f.readlines())
	return observed

observed= read_in(training_file)

# data structures required for algo 
L = len(observed)
interval = [0]*L
PI={}
count=0

with open(result_training) as f:
	while 1:
		s=f.readline()
		if s=='':
			break
		a,b=map(int,s.split())
		interval[a]+=1
		interval[b+1]-=1

for x in range(1,L):
	interval[x]+=interval[x-1]

hidden = list(observed)
for x in range(1,L):
	hidden[x]=hidden[x]+('+' if interval[x]>0 else '-')
	if interval[x]==0:
		PI[hidden[x]]=PI.get(hidden[x],0)+1.0
		count+=1
	if hidden[x] not in PI:
		PI[hidden[x]]=OWNZERO

for x in PI:
	PI[x]=PI[x]*1.0/count


freq={}
transition={}
emission={}
for x in range(1,L-1):
	freq[ hidden[x] ]=freq.get(hidden[x],0)+1
	if hidden[x] not in transition:
		transition[hidden[x]]={}
	if hidden[x] not in emission:
		emission[hidden[x]]={}
	transition[hidden[x]][hidden[x+1]]= transition[hidden[x]].get(hidden[x+1],0)+1

for x in emission:
	for y in PI:
		emission[x][y[0]]=1.0 - OWNZERO if y[0]==x[0] else OWNZERO
		if y not in transition[x]:
			transition[x][y]=OWNZERO
	# 	print emission[x][y[0]],
	# print
	
for x in transition:
	# print transition[x]
	for y in transition[x]:
		transition[x][y]*=1.0/freq[x]
	# 	print transition[x][y],
	# print

# print PI

# VITERBI IMPLEMENTATION.
def Viterbi(obs_seq):
	length = len(obs_seq)
	V= [ {x:float("-inf") for x in PI }  for i in range(length)  ]
	back= [ {x:x for x in PI }  for i in range(length)  ]
	for j in PI:
		V[1][j]=log(emission[j][obs_seq[1]]) + log(PI[j]) 

	for i in range(2,length):
		for j in PI:
			best,trace=float("-inf"),''
			for k in PI:
				c = V[i-1][k] + log(transition[k][j]) + log(emission[j][obs_seq[i]])
				if best<c:
					best=c
					trace=k
			V[i][j]=best
			back[i][j]=trace

	output=[' ' for i in range(length)]

	print( "log-liklihood:",max(V[length-1].items(), key=operator.itemgetter(1))[1])
	output[length-1]=max(V[length-1].items(), key=operator.itemgetter(1))[0]
	for x in range(length-2,0,-1):
		j=back[x+1][output[x+1]]
		output[x]=j
	return output

def write_out(seq):
	cpgs=[]
	i=1
	length=len(seq)
	while i<length:
		if seq[i][1]=='+':
			j=i+1
			while(j<length):
				if(seq[j][1]=='-'):
					break
				j+=1
			cpgs.append([i,j-1])
			i=j
		i+=1
	count = len(cpgs)
	return count


start_time=time.time()