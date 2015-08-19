#!/usr/bin/python
import os
import shutil

nproc = 25

dirs = []

l = []

if (nproc > 9):
	for i in range(0,10):
		name = 'map-0' + str(i)
		dirs.append(name)
	for i in range(10,nproc):
		name = 'map-'  + str(i)
		dirs.append(name)
else:
	for i in range(0,nproc):
		name = 'map-'  + str(i)
		dirs.append(name)

currentpath = os.getcwd()
for i in range(0,nproc):
	workingd = currentpath + '/' + dirs[i]
	os.chdir(workingd)
	os.system('qsub submit_cluster.pbs')
