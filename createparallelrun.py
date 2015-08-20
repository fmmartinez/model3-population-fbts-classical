#!/usr/bin/python
import os
import shutil

nproc = 10

dirs = []

l = []
m = []

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

#generate folders
for i in range(0,nproc):
	if not(os.path.exists(dirs[i])):
			os.mkdir(dirs[i])

#original structure of md.in file in map00:
l.append('Np\tDELTA\tNOSC\tOME_MAX\t\r\n')
l.append('1\t1d0\t20\t50\r\n')
l.append('NMCS\tNMDS\tseed\tDT\tLAMDA_D\r\n')
l.append('5000\t4000\t5\t5d-5\t10d0\r\n')
l.append('Eg\tEb\tEd\tmu\tE0\tbeta\r\n')
l.append('0\t240\t240\t1\t0.90\t0.24\r\n')
l.append('TIME_J\tTAU_J\tOMEGA_J\r\n')
l.append('0.15d0\t0.045d0\t260d0\r\n')
l.append('BATH(0:B EQ 1:T EQ)\tINIT\r\n')
l.append('0\t3\r\n')
l.append('basispercenter\tmapinG\tmapinB\tmapinD\tcount\r\n')
l.append('1\t1\t1\t1\t0\r\n')

#creating md.in files
for i in range(0,nproc):
	mdfile = open('./' + dirs[i] + '/md.in','w')
	
	seed = i*9 + 5

	l[3] = '500000\t4000\t' + str(seed) + '\t5d-5\t10d0\r\n'

	for j in range(0,12):
		mdfile.write(l[j])

	mdfile.close()

#original structure of pbs fie
m.append('#!/bin/bash -l\r\n')
m.append('#PBS -S /bin/bash\r\n')
m.append('#PBS -N pbme\r\n')
m.append('#PBS -l walltime=100:00:00\r\n')
m.append('cd $PBS_O_WORKDIR\r\n')
m.append('time ./a.out < md.in\r\n')

#creating pbs files
for i in range(0,nproc):
	pbsfile = open('./' + dirs[i] + '/submit_cluster.pbs','w')
	
	m[2] = '#PBS -N fbc-' + str(i) + '\r\n' 

	for j in range(0,6):
		pbsfile.write(m[j])

	pbsfile.close()


#copy executables
for i in range(0,nproc):
	shutil.copy2('a.out',dirs[i])
