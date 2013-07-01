# Produces a filtered point process MER using the rate data given from the neural mass simulation given in BrainSim.py. Requires 2 arguments; number of neurons and STN rate data file.
#
# Kristian Weegink: uqkweegi@uq.edu.au
# 07/2013

import sys
import numpy as np
from matplotlib import pylab
import random

if (len(sys.argv)<2 or int(sys.argv[1])<1):
	N=10000
else:
	N=int(sys.argv[1])

if (len(sys.argv)<3):
	STNrate='h:\data\BrainSimSTN'
else:
	STNrate=sys.argv[2]
	
data = np.loadtxt(STNrate,delimiter=' ')
print "Rate data loaded"

STNdata=[]
tick=[]

for n in data:
	STNdata.append(n[1])
	tick.append(n[0])

Ratetime=pylab.cumsum(tick)
	
maxrate=1. /0.009
dt=1. /24000

times=[]

for n in range(int(Ratetime[-1]/dt)):
	times.append(dt*n)

It=np.loadtxt('h:\data\matlab\\apcurrent24k.dat',delimiter=',')
print 'Current loaded'

It=np.multiply(np.true_divide(It,It.min()),250e-9)          #normalize
currentLength=len(It)-1
epsilon=8.85e-12                      #Permitivity of free space
rho=10.**5 * 10.**6        	   	  		#density of neurons in STN m^-3
r=np.power(np.multiply(3./4*N/(np.pi*rho),np.array([random.uniform(0,1) for _ in range(N)])),1./3)   #create a power law distribution of neuron radii

r.sort()

R3=0.96e3
C3=2.22e-6
C2=9.38e-9
C3=1.56e-6
C2=9.38e-9
R4=100.e6
R2N=np.multiply(1./(4*np.pi*epsilon),r[-1:])
R1=2100.;
t_impulse=np.array([dt*n for n in range(100)])

print 'initialization complete'

for neuron in range(N):
	