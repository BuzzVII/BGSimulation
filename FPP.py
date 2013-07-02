# Produces a filtered point process MER using the rate data given from the neural mass simulation given in BrainSim.py. Requires 2 arguments; number of neurons and STN rate data file.
#
# Kristian Weegink: uqkweegi@uq.edu.au
# 07/2013

import sys
import numpy as np
from scipy import signal
from matplotlib import pylab
import random

print '\n'*100

if (len(sys.argv)<2 or int(sys.argv[1])<1):
	N=10000
else:
	N=int(sys.argv[1])

if (len(sys.argv)<3):
	STNrate='C:\\Users\\uqkweegi\\Dropbox\\Code\\Matlab\\PPpaper\\STN'#h:\data\BrainSimSTN'
else:
	STNrate=sys.argv[2]
	
Ctxrate='C:\\Users\\uqkweegi\\Dropbox\\Code\\Matlab\\PPpaper\\ctx'

data = np.loadtxt(STNrate,delimiter=' ')
datactx=np.loadtxt(Ctxrate,delimiter=' ')
print "Rate data loaded"

STNdata=[]
Ctxdata=[]
tick=[]

for n in data:
	STNdata.append(n[1])
	tick.append(n[0])

for n in datactx:
	Ctxdata.append(n[1])
	
Ratetime=pylab.cumsum(tick)
BGdt=tick[1]
	
maxrate=1. /0.009
dt=1. /24000

times=[]

for n in range(int(Ratetime[-1]/dt)):
	times.append(dt*n)

It=np.loadtxt('C:\\Users\\uqkweegi\\Dropbox\\Code\\Matlab\\PPpaper\\apcurrent24k.dat',delimiter=',') #h:\data\matlab\\apcurrent24k.dat',delimiter=',')
print 'Current loaded'

It=np.multiply(np.true_divide(It,It.min()),250e-9)          #normalize
currentLength=len(It)
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
R2N=np.multiply(1./(4*np.pi*epsilon),r)
R1=2100.;
t_impulse=np.array([dt*n for n in range(100)])

print 'initialization complete'

Vt=pylab.zeros(len(times))

for neuron in range(N):
	R2=R2N[neuron]
	ppwave=pylab.zeros(len(times))
	absoluteTimes=np.random.exponential(1./(maxrate*STNdata[0]),1)
	while absoluteTimes[-1] < times[-1]-currentLength*dt:
		wave_start=int(absoluteTimes[-1]/dt)
		wave_end=wave_start+currentLength
		if wave_end > len(times):
			break
		ppwave[wave_start:wave_end]=np.add(ppwave[wave_start:wave_end],It)
		isi=np.random.exponential(1./(maxrate*STNdata[int(absoluteTimes[-1]/BGdt)]),1)
		absoluteTimes=np.append(absoluteTimes,[absoluteTimes[-1]+isi])
##############################################
	extracellular_impulse_response=np.multiply(np.multiply(np.exp(np.multiply(t_impulse,-((C2*R1*R2 + C2*R1*R3 + C2*R1*R4 - C3*R1*R3 + C3*R2*R3 + C3*R3*R4))/(2*C2*C3*R1*R3*(R2 + R4)))),(np.add(np.cosh(np.multiply(t_impulse,(C2**2*R1**2*R2**2 + 2*C2**2*R1**2*R2*R3 + 2*C2**2*R1**2*R2*R4 + C2**2*R1**2*R3**2 + 2*C2**2*R1**2*R3*R4 + C2**2*R1**2*R4**2 + 2*C2*C3*R1**2*R2*R3 - 2*C2*C3*R1**2*R3**2 + 2*C2*C3*R1**2*R3*R4 - 2*C2*C3*R1*R2**2*R3 - 2*C2*C3*R1*R2*R3**2 - 4*C2*C3*R1*R2*R3*R4 - 2*C2*C3*R1*R3**2*R4 - 2*C2*C3*R1*R3*R4**2 + C3**2*R1**2*R3**2 - 2*C3**2*R1*R2*R3**2 - 2*C3**2*R1*R3**2*R4 + C3**2*R2**2*R3**2 + 2*C3**2*R2*R3**2*R4 + C3**2*R3**2*R4**2)**(1/2)/(2*C2*C3*R1*R3*(R2 + R4)))),np.divide(np.sinh(np.multiply(t_impulse,(C2**2*R1**2*R2**2 + 2*C2**2*R1**2*R2*R3 + 2*C2**2*R1**2*R2*R4 + C2**2*R1**2*R3**2 + 2*C2**2*R1**2*R3*R4 + C2**2*R1**2*R4**2 + 2*C2*C3*R1**2*R2*R3 - 2*C2*C3*R1**2*R3**2 + 2*C2*C3*R1**2*R3*R4 - 2*C2*C3*R1*R2**2*R3 - 2*C2*C3*R1*R2*R3**2 - 4*C2*C3*R1*R2*R3*R4 - 2*C2*C3*R1*R3**2*R4 - 2*C2*C3*R1*R3*R4**2 + C3**2*R1**2*R3**2 - 2*C3**2*R1*R2*R3**2 - 2*C3**2*R1*R3**2*R4 + C3**2*R2**2*R3**2 + 2*C3**2*R2*R3**2*R4 + C3**2*R3**2*R4**2)**(1/2)/(2*C2*C3*R1*R3*(R2 + R4))))*(C2*R1*R2 - C2*R1*R3 + C2*R1*R4 + C3*R1*R3 - C3*R2*R3 - C3*R3*R4),(C2**2*R1**2*R2**2 + 2*C2**2*R1**2*R2*R3 + 2*C2**2*R1**2*R2*R4 + C2**2*R1**2*R3**2 + 2*C2**2*R1**2*R3*R4 + C2**2*R1**2*R4**2 + 2*C2*C3*R1**2*R2*R3 - 2*C2*C3*R1**2*R3**2 + 2*C2*C3*R1**2*R3*R4 - 2*C2*C3*R1*R2**2*R3 - 2*C2*C3*R1*R2*R3**2 - 4*C2*C3*R1*R2*R3*R4 - 2*C2*C3*R1*R3**2*R4 - 2*C2*C3*R1*R3*R4**2 + C3**2*R1**2*R3**2 - 2*C3**2*R1*R2*R3**2 - 2*C3**2*R1*R3**2*R4 + C3**2*R2**2*R3**2 + 2*C3**2*R2*R3**2*R4 + C3**2*R3**2*R4**2)**(1/2))))),-R4/(C2*(R2 + R4)));
	electrode_ppwave=np.convolve(ppwave,extracellular_impulse_response,'same');
########################################	
	Vt=np.add(Vt,electrode_ppwave)
	if np.mod(neuron,1000)==999:
		print(str(neuron+1)+" neurons calculated")
		
print 'neuron contribution to MER complete'
	
Vt=np.subtract(Vt,np.mean(Vt))

flow=10000.;
fhigh=500.;

b,a=signal.butter(9,flow/24000,'low');
Vt=signal.lfilter(b, a, Vt);
b,a=signal.butter(1,fhigh/24000,'high');
Vt=signal.lfilter(b, a, Vt);

volts=pylab.plot(times,Vt)
stnrate=pylab.plot(Ratetime,np.multiply(STNdata,200))
pylab.show()