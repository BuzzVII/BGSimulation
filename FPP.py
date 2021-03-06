#!/usr/bin/python 
#
# Produces a filtered point process MER using the rate data given from the neural mass simulation given in BrainSim.py. Requires 2 arguments; number of neurons and STN rate data file.
#
# Kristian Weegink: uqkweegi@uq.edu.au
# 07/2013

import sys
import numpy as np
from scipy import signal
from matplotlib import pylab
import random
import Logger
import math

def FPP(log=Logger.logger(0), N = 10000, dt = 1./24000, distributionParameter = [30], plotAll = True, efield = False):
    
    #check if rate file or rate is present
    if len(distributionParameter)  ==  1:
        try: 
            data = np.loadtxt(distributionParameter[0],delimiter = ' ')
            log.info("Rate data loaded")
            BGsim = True
            STNdata = []
            tick = []
            for n in data:
                STNdata.append(n[1])
                tick.append(n[0])
            Ratetime = pylab.cumsum(tick)
            BGdt = tick[1]
            timeSteps = int(Ratetime[-1]/dt)  
        except:
            float(distributionParameter[0])
            Ratetime = 1.
            timeSteps = int(Ratetime/dt)
            BGsim = False
    else:
        Ratetime = 1.
        BGsim = False
        timeSteps = int(Ratetime/dt)
        
    maxrate = 1./0.009
    times = []
    for n in range(timeSteps):
        times.append(dt*n)

    # check for current file, if none present use impules
    try:
        It = np.loadtxt('C:\\Users\\Kristian\\Dropbox\\phd\\Data\\apcurrent24k.dat',delimiter = ',')

        #/home/uqkweegi/Documents/Data/apcurrent24k.dat',delimiter = ',')
    except:
        log.error('no current file present')
        It = np.array(1)

    log.info('Current loaded')
    It = np.multiply(np.true_divide(It,It.min()),250e-9)          #normalize
    currentLength = len(It)
    
    #calculate extracellular effects
    epsilon = 8.85e-12                                #Permitivity of free space
    rho = 10.**5 * 10.**6                     #density of neurons in STN m^-3
    r = np.power(np.multiply(3./4*N/(np.pi*rho),np.array([random.uniform(0,1) for _ in range(N)])),1./3)   #create a power law distribution of neuron radii
    r.sort()
    if efield:
        rijk = [[random.uniform(0,1)-0.5 for _ in range(N)],[random.uniform(0,1)-0.5 for _ in range(N)],[random.uniform(0,1)-0.5 for _ in range(N)]] #create vector direction of field 
        #if plotAll:
        #    vi = pylab.plot(rijk[0])
        #    vj = pylab.plot(rijk[1])
        #    vk = pylab.plot(rijk[2])
        #    pylab.show()
    R3 = 0.96e3
    C3 = 2.22e-6
    C2 = 9.38e-9
    C3 = 1.56e-6
    C2 = 9.38e-9
    R4 = 100.e6
    R2N = np.multiply(1./(4*np.pi*epsilon),r)
    R1 = 2100.;
    t_impulse = np.array([dt*n for n in range(100)])

    log.info('initialization complete')

    Vt = pylab.zeros(len(times))
    Vi = Vt
    Vj = Vt
    Vk = Vt

    # start simulation
    #-------------------------------------------------------------------------------#
    for neuron in range(N):
        R2 = R2N[neuron]
        ppwave = pylab.zeros(len(times))
        if BGsim:
            absoluteTimes = np.random.exponential(1./(maxrate*STNdata[0]),1)
        else:
            if len(distributionParameter)  ==  1:
                absoluteTimes = np.random.exponential(1./(distributionParameter[0]),1)
            else:
                absoluteTimes = [random.weibullvariate(distributionParameter[0],distributionParameter[1])]
        while absoluteTimes[-1] < times[-1]-currentLength*dt:
            wave_start = int(absoluteTimes[-1]/dt)
            wave_end = wave_start+currentLength
            if wave_end > len(times):
                break
            ppwave[wave_start:wave_end] = np.add(ppwave[wave_start:wave_end],It)
            if BGsim:
                isi = np.random.exponential(1./(maxrate*STNdata[int(absoluteTimes[-1]/BGdt)]),1)
            else:
                if len(distributionParameter)  ==  1:
                    isi = np.random.exponential(1./(distributionParameter[0]),1)
                else:
                    isi = random.weibullvariate(distributionParameter[0],distributionParameter[1])
            absoluteTimes = np.append(absoluteTimes,[absoluteTimes[-1]+isi])
        # calculate neuron contribution
        #------------------------------------------------------------------------------#
        extracellular_impulse_response = np.multiply(np.multiply(np.exp(np.multiply(t_impulse,-20*17*((C2*R1*R2 + C2*R1*R3 + C2*R1*R4 - C3*R1*R3 + C3*R2*R3 + C3*R3*R4))/(2*C2*C3*R1*R3*(R2 + R4)))),(np.add(np.cosh(np.multiply(t_impulse,(C2**2*R1**2*R2**2 + 2*C2**2*R1**2*R2*R3 + 2*C2**2*R1**2*R2*R4 + C2**2*R1**2*R3**2 + 2*C2**2*R1**2*R3*R4 + C2**2*R1**2*R4**2 + 2*C2*C3*R1**2*R2*R3 - 2*C2*C3*R1**2*R3**2 + 2*C2*C3*R1**2*R3*R4 - 2*C2*C3*R1*R2**2*R3 - 2*C2*C3*R1*R2*R3**2 - 4*C2*C3*R1*R2*R3*R4 - 2*C2*C3*R1*R3**2*R4 - 2*C2*C3*R1*R3*R4**2 + C3**2*R1**2*R3**2 - 2*C3**2*R1*R2*R3**2 - 2*C3**2*R1*R3**2*R4 + C3**2*R2**2*R3**2 + 2*C3**2*R2*R3**2*R4 + C3**2*R3**2*R4**2)**(1/2)/(2*C2*C3*R1*R3*(R2 + R4)))),np.divide(np.sinh(np.multiply(t_impulse,(C2**2*R1**2*R2**2 + 2*C2**2*R1**2*R2*R3 + 2*C2**2*R1**2*R2*R4 + C2**2*R1**2*R3**2 + 2*C2**2*R1**2*R3*R4 + C2**2*R1**2*R4**2 + 2*C2*C3*R1**2*R2*R3 - 2*C2*C3*R1**2*R3**2 + 2*C2*C3*R1**2*R3*R4 - 2*C2*C3*R1*R2**2*R3 - 2*C2*C3*R1*R2*R3**2 - 4*C2*C3*R1*R2*R3*R4 - 2*C2*C3*R1*R3**2*R4 - 2*C2*C3*R1*R3*R4**2 + C3**2*R1**2*R3**2 - 2*C3**2*R1*R2*R3**2 - 2*C3**2*R1*R3**2*R4 + C3**2*R2**2*R3**2 + 2*C3**2*R2*R3**2*R4 + C3**2*R3**2*R4**2)**(1/2)/(2*C2*C3*R1*R3*(R2 + R4))))*(C2*R1*R2 - C2*R1*R3 + C2*R1*R4 + C3*R1*R3 - C3*R2*R3 - C3*R3*R4),(C2**2*R1**2*R2**2 + 2*C2**2*R1**2*R2*R3 + 2*C2**2*R1**2*R2*R4 + C2**2*R1**2*R3**2 + 2*C2**2*R1**2*R3*R4 + C2**2*R1**2*R4**2 + 2*C2*C3*R1**2*R2*R3 - 2*C2*C3*R1**2*R3**2 + 2*C2*C3*R1**2*R3*R4 - 2*C2*C3*R1*R2**2*R3 - 2*C2*C3*R1*R2*R3**2 - 4*C2*C3*R1*R2*R3*R4 - 2*C2*C3*R1*R3**2*R4 - 2*C2*C3*R1*R3*R4**2 + C3**2*R1**2*R3**2 - 2*C3**2*R1*R2*R3**2 - 2*C3**2*R1*R3**2*R4 + C3**2*R2**2*R3**2 + 2*C3**2*R2*R3**2*R4 + C3**2*R3**2*R4**2)**(1/2))))),-R4/(C2*(R2 + R4)));
        electrode_ppwave = np.convolve(ppwave,extracellular_impulse_response,'same');
        if efield:  #add fields
            amp = 1/np.sqrt((np.square(rijk[0][neuron])+np.square(rijk[1][neuron])+np.square(rijk[2][neuron])))
            rijk[0][neuron] = rijk[0][neuron]*amp
            rijk[1][neuron] = rijk[1][neuron]*amp
            rijk[2][neuron] = rijk[2][neuron]*amp
            Vi = np.add(Vi,np.multiply(electrode_ppwave,rijk[0][neuron]))
            Vj = np.add(Vj,np.multiply(electrode_ppwave,rijk[1][neuron]))
            Vk = np.add(Vk,np.multiply(electrode_ppwave,rijk[2][neuron]))
        else:       #add scalar
            Vt = np.add(Vt,electrode_ppwave)
        if np.mod(neuron,1000) == 999:
            log.info(str(neuron+1)+" neurons calculated")
    #------------------------------------------------------------------------------#        
    # end simulation
    
    log.info('neuron contribution to MER complete')
    
    #remove bias
    if efield:
        Vt = np.sqrt(np.add(np.square(Vi),np.square(Vj),np.square(Vk)))   
    Vt = np.subtract(Vt,np.mean(Vt))

    #apply hardware filters
    flow = 5500*2.
    fhigh = 500.
    b,a = signal.butter(18,flow*dt,'low')
    Vt = signal.lfilter(b, a, Vt)
    b,a = signal.butter(1,fhigh*dt,'high')
    Vt = signal.lfilter(b, a, Vt)

    #produce plots
    if plotAll:
        volts = pylab.plot(times,Vt)
        if BGsim:
            stnrate = pylab.plot(Ratetime,np.multiply(STNdata,200))
        pylab.show()
        nfft=2**int(math.log(len(Vt),2))+1
        sr = 1/dt
        Pxi,freqs=pylab.psd(x=Vt,Fs=sr,NFFT=nfft/10,window=pylab.window_none, noverlap=100)
        pylab.show()
        return freqs, Pxi
        psd = pylab.loglog(freqs, Pxi)
        pylab.show()
    return Vt, times
    
def main():

    print '\n'*100      #clear screen
    
    log = Logger.logger(loglevel=1)

    if (len(sys.argv)<2 or int(sys.argv[1])<1):
        N = 10000
    else:
        N = int(sys.argv[1])

    if (len(sys.argv)<4):
        STNrate = ['C:\\Users\\Kristian\\Dropbox\\phd\\Data\\STN']
    else:
        STNrate = [float(sys.argv[3])]

    #if (len(sys.argv)<3):
    #    field = sys.argv[2]=='True'
    #else:
    field = False
 
    Vt,times = FPP(log,N, distributionParameter = STNrate, efield = field, plotAll = True)



if (__name__  ==  '__main__'):
    main()
