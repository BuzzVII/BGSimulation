# wrapper for an object that opens a patient file finds the power specturm
# then runs a simulation given C and rate and finds the returns the residues
# of the difference of the PSDs
#
# Kristian Weegink, Nov 2013

import numpy as np
from scipy import signal
from matplotlib import pylab
import FPP
import BrainSim as BS
import math
import scipy.io.wavfile as wav
import Logger as lg

class Patient:
    
    # intialize by finding the patient PSD
    def __init__(self,patient_file,logger=lg.logger(0)):
        self.log = logger
        self.sr,self.Vt = wav.read(patient_file)
        self.log.info('patient loaded')
        self.nfft=2**int(math.log(len(self.Vt),2))+1
        self.Pxx,self.freqs=pylab.psd(x=self.Vt,Fs=self.sr,NFFT=self.nfft,window=pylab.window_none, noverlap=1)   
        self.log.info('patient initialized')

    # find the L2 difference between the FPP simulation and patient PSDs
    def fpp(self,rate,C):
        self.log.info('C = ' + str(C) + '; rate = ' + str(rate) )
        Vt,t = FPP.FPP(log=lg.logger(0), N=10000, dt=1./24000, distributionParameter=[C,rate], plotAll=False,efield=False)
        Pxi,freqs = pylab.psd(x=Vt,Fs=self.sr,NFFT=self.nfft,window=pylab.window_none, noverlap=1)
        self.log.info('L2  = ' + str(np.sqrt((np.square(np.subtract(self.Pxx[370:3300],10000*Pxi[370:3300])).sum()))))
        #pylab.show()
        #pylab.loglog(10000*Pxi[37:3300])
        #pylab.loglog(self.Pxx[37:3300])
        #pylab.show()
        #pylab.plot(self.freqs)
        #pylab.show()
        return np.subtract(self.Pxx[37:3300],10000*Pxi[37:3300])
