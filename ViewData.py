# Script to view output from BrainSim.py neural mass simulator. Based on pprateshape.m in FPP project
#
# Kristian Weegink: uqkweegi@uq.edu.au
# 07/2013 

import numpy as np
from matplotlib import pylab

data = np.loadtxt('C:\\Users\\Kristian\\Dropbox\\phd\\Data\\STN',delimiter=' ')

STNdata=[]
tick=[]

for n in data:
	STNdata.append(n[1])
	tick.append(n[0])

data = np.loadtxt('C:\\Users\\Kristian\\Dropbox\\phd\\Data\\GPe',delimiter=' ')

GPedata=[]

for n in data:
	GPedata.append(n[1])
	
time=pylab.cumsum(tick)
time=time*1000
	
figure1=pylab.plot(time, STNdata,label='STN')
figure2=pylab.plot(time, GPedata,label='GPe')

pylab.legend(loc=4)
pylab.xlabel('time (ms)')
pylab.ylabel('fraction of max rate')

pylab.show()
